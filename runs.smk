#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download data from SRA
"""
localrules: sra_download

rule sra_download:
    output:
        "runs/{run_acc}/unmapped.bam"
    log:
        "runs/{run_acc}/fasterq-dump.log",
        "runs/{run_acc}/fastq_to_sam.log",
        "runs/{run_acc}/mark_adapters.log",
        "runs/{run_acc}/mark_adapters.metrics"
    params: 
        keyfile = config['dbgap_key'] if 'dbgap_key' in config else 'nonexistent',
        tmpdir = config['tmpdir'],
        rundir = "runs/{run_acc}"
    conda:
        "envs/utils.yaml"    
    threads: min(20, snakemake.utils.available_cpu_count())
    shell:
        '''
        mkdir -p {params.rundir}
        
        #--- Download using fasterq-dump
        if [[ ! -e "{params.keyfile}" ]]; then
            fasterq-dump\
              -e {threads}\
              -t {params.tmpdir}\
              -O {params.rundir}\
              {wildcards.run_acc}\
            &> {log[0]}
        else
            CDIR=$PWD
            NCBIDIR=$(vdb-config --import {params.keyfile} | tail -n+2 | cut -d"'" -f2)
            cd $NCBIDIR
            fasterq-dump\
              -e {threads}\
              -t {params.tmpdir}\
              {wildcards.run_acc}\
            &> $CDIR/{log[0]}
            cd $CDIR
            mv $NCBIDIR/{wildcards.run_acc}* {params.rundir}/
        fi
        
        #--- FASTQ to uBAM with marked adapters
        if [ -e {params.rundir}/{wildcards.run_acc}_1.fastq ]; then
            picard FastqToSam\
              TMP_DIR={params.tmpdir}\
              F1={params.rundir}/{wildcards.run_acc}_1.fastq\
              F2={params.rundir}/{wildcards.run_acc}_2.fastq\
              O=/dev/stdout\
              SM={wildcards.run_acc}\
              RG={wildcards.run_acc}\
            2> {log[1]} |\
            picard MarkIlluminaAdapters\
              I=/dev/stdin O={output[0]} M={log[3]}\
              COMPRESSION_LEVEL=5 TMP_DIR={params.tmpdir}\
            2> {log[2]}
        else
            picard FastqToSam\
              TMP_DIR={params.tmpdir}\
              F1={params.rundir}/{wildcards.run_acc}.fastq\
              O=/dev/stdout\
              SM={wildcards.run_acc}\
              RG={wildcards.run_acc}\
            2> {log[1]} |\
            picard MarkIlluminaAdapters\
              I=/dev/stdin O={output[0]} M={log[3]}\
              COMPRESSION_LEVEL=5 TMP_DIR={params.tmpdir}\
            2> {log[2]}
        fi
        
        chmod 600 {output[0]}
        echo "rm -f {params.rundir}/{wildcards.run_acc}*.fastq"
        '''
    


"""
Revert mapped data to unmapped
"""

# Default SAM attributes cleared by RevertSam
attr_revertsam = ['NM', 'UQ', 'PG', 'MD', 'MQ', 'SA', 'MC', 'AS']
# SAM attributes output by STAR
attr_star = ['NH', 'HI', 'NM', 'MD', 'AS', 'nM', 'jM', 'jI', 'XS', 'uT']
# Additional attributes to clear
ALN_ATTRIBUTES = list(set(attr_star) - set(attr_revertsam))


rule revert_and_mark_adapters:
    input:
        "runs/{run_acc}/original.bam"
    output:
        "runs/{run_acc}/unmapped.bam"
    log:
        "runs/{run_acc}/revert_bam.log",
        "runs/{run_acc}/mark_adapters.log",
        "runs/{run_acc}/mark_adapters.metrics"
    conda:
        "envs/utils.yaml"
    params:
        attr_to_clear = expand("ATTRIBUTE_TO_CLEAR={a}", a=ALN_ATTRIBUTES)
    shell:
       '''
        picard RevertSam\
          I={input[0]} O=/dev/stdout\
          SANITIZE=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT\
          TMP_DIR={config[tmpdir]} {params.attr_to_clear}\
        2> {log[0]} |\
        picard MarkIlluminaAdapters\
          I=/dev/stdin O={output[0]} M={log[2]}\
          COMPRESSION_LEVEL=5 TMP_DIR={config[tmpdir]}\
        2> {log[1]}
        
        chmod 600 {output[0]}
        '''


"""
Align with bowtie2 allowing for multiple alignments
"""
rule bowtie2_multi:
    input:
        "runs/{run_acc}/unmapped.bam",
        ancient(expand(config['indexes']['bowtie2'] + ".{i}.bt2", i = range(1, 5))),
        ancient(expand(config['indexes']['bowtie2'] + ".rev.{i}.bt2", i = range(1, 3)))        
    output:
        "runs/{run_acc}/aligned.bt2.bam"
    log:
        "runs/{run_acc}/aligned.bt2.log"
    conda:
        "envs/bowtie2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.run_acc}.XXXXXX)
        picard SamToFastq\
          I={input[0]} F=$tdir/R1.fq F2=$tdir/R2.fq\
          CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 NON_PF=true TMP_DIR=$tdir

        bowtie2 -p {threads} -X 1200\
          -x {config[indexes][bowtie2]}\
          -1 $tdir/R1.fq -2 $tdir/R2.fq\
          --rg-id {wildcards.run_acc} --rg PL:ILLUMINA \
          -k 100 --very-sensitive-local --score-min L,0,1.6\
        2> {log[0]} | samtools view -b > {output[0]}

        chmod 600 {output[0]}
        rm -rf $tdir
        '''


localrules: locate_sra
rule locate_sra:
    input:
        config['sra_run_table']
    output:
        '%s.URL.txt' % config['sra_run_table']
    params: 
        keyfile = config['dbgap_key'] if 'dbgap_key' in config else None
    run:
        out_header = ['accession', 'filename', 'url', 'size', 'md5']
        out_table = []
        
        # Construct query
        accessions = [d[config['colname_RUN']] for d in METADATA]
        request_params = [('acc', a) for a in accessions]
        request_params.append(('location', config['sra_location']))
        request_params.append(('type', 'sra'))
        
        # Authorization file
        if params.keyfile is not None and os.path.exists(params.keyfile):
            files = {'ngc': (params.keyfile, open(params.keyfile, 'rb'))}
        else:
            files = None

        # Send request to SRA data locator
        locator_url = 'https://www.ncbi.nlm.nih.gov/Traces/sdl/1/retrieve'
        response = requests.post(locator_url, params=request_params, files=files)
        for r in response.json():
            if r['status'] == 200:
                for f in r['files']:
                    out_table.append([
                        r['accession'], f['name'], f['link'], f['size'], f['md5'], 
                    ])
            else:
                print('%s\tERROR[%s]\t%s' % (r['accession'], r['status'], r['message']))
        # Print output table
        if out_table:
            with open(output[0], 'w') as outh:
                print('\t'.join(out_header), file=outh)
                print('\n'.join('\t'.join(r) for r in out_table), file=outh)
