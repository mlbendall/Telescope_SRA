#! /usr/bin/env python
# -*- coding: utf-8 -*-

localrules: combine_bams
rule combine_bams:
    input:
        bams = lambda wildcards: expand("runs/{f}/{p}.bam", 
                                        p=wildcards.bampre,
                                        f=SAMPLE_RUN[wildcards.sampid])
    output:
        "samples/{sampid}/{bampre}.bam"
    params:
        nfiles = lambda wildcards: len(SAMPLE_RUN[wildcards.sampid])    
    conda:
        "envs/utils.yaml"
    shell:
        '''
        mkdir -p $(dirname {output[0]})
        if [[ "1" -eq "{params.nfiles}" ]]; then
            ln -s ../../{input.bams[0]} {output[0]}
        else
            tmphead=$(dirname {output[0]})/header.sam
            scripts/merge_sam_headers.py {input.bams} > $tmphead
            samtools cat -h $tmphead -o {output[0]} {input.bams}
            rm -f $tmphead
        fi
        chmod 600 {output[0]}
        '''


localrules: combine_bt2_logs
rule combine_bt2_logs:
    input:
        logs = lambda wildcards: expand("runs/{f}/aligned.bt2.log", 
                                        f=SAMPLE_RUN[wildcards.sampid])
    output:
        "samples/{sampid}/aligned.bt2.log"
    shell:
        'scripts/combine_bowtie2_logs.py {input.logs} > {output[0]}'


rule kallisto:
    input:
        "samples/{sampid}/unmapped.bam",
        config['indexes']['kallisto']  
    output:
        "samples/{sampid}/kallisto.abundance.h5",
        "samples/{sampid}/kallisto.abundance.tsv",
        "samples/{sampid}/kallisto.run_info.json" 
    log:
        "samples/{sampid}/kallisto.log"
    conda:
        "envs/kallisto.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)

        picard SamToFastq\
          I={input[0]} F=$tdir/R1.fq F2=$tdir/R2.fq\
          CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=N NON_PF=true TMP_DIR=$tdir
        
        HDF5_USE_FILE_LOCKING=FALSE kallisto quant\
          -t {threads} -b 100 -i {input[1]} -o $tdir\
          $tdir/R1.fq $tdir/R2.fq 2>&1 | tee {log[0]}

        mv $tdir/abundance.h5 {output[0]}
        mv $tdir/abundance.tsv {output[1]}
        mv $tdir/run_info.json {output[2]}
        
        rm -rf $tdir
        '''


rule telescope:
    input:
        aln = "samples/{sampid}/aligned.bt2.bam",
        ann = config['annotations']['retro']
    output:
        "samples/{sampid}/telescope.report.tsv",
        "samples/{sampid}/telescope.updated.bam"
    log:
        "samples/{sampid}/telescope.log"
    conda:
        "envs/telescope.yaml"
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.sampid}.XXXXXX)
        
        telescope assign\
          --exp_tag inform --theta_prior 200000 --max_iter 200 --updated_sam\
          --outdir $tdir {input[0]} {input[1]} 2>&1 | tee {log[0]}
        
        mv $tdir/inform-telescope_report.tsv {output[0]}
        mv $tdir/inform-updated.bam {output[1]}
        chmod 600 {output[1]}
        
        rm -rf $tdir
        '''


localrules: sample_complete
rule sample_complete:
    input:
        rules.kallisto.output,
        rules.telescope.output
    output:
        touch("samples/{sampid}/completed.txt")

