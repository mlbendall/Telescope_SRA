#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Download data from GDC
"""
localrules: gdc_download

rule gdc_download:
    output:
        "files/{uuid}/original.bam"
    params:
        md5sum = lambda wildcards: FILE_MD5[wildcards.uuid]
    shell:
        '''
        mkdir -p $(dirname {output[0]})
        curl --header "X-Auth-Token: $(<{config[gdc_token_file]})"\
          https://api.gdc.cancer.gov/data/{wildcards.uuid} > {output[0]}
        chmod 600 {output[0]}
        echo {params.md5sum} {output[0]} | md5sum -c -
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
        "files/{uuid}/original.bam"
    output:
        "files/{uuid}/unmapped.bam"
    log:
        "files/{uuid}/revert_bam.log",
        "files/{uuid}/mark_adapters.log",
        "files/{uuid}/mark_adapters.metrics"
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
        "files/{uuid}/unmapped.bam",
        ancient(expand(config['indexes']['bowtie2'] + ".{i}.bt2", i = range(1, 5))),
        ancient(expand(config['indexes']['bowtie2'] + ".rev.{i}.bt2", i = range(1, 3)))        
    output:
        "files/{uuid}/aligned.bt2.bam"
    log:
        "files/{uuid}/aligned.bt2.log"
    conda:
        "envs/bowtie2.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
        tdir=$(mktemp -d {config[tmpdir]}/{rule}.{wildcards.uuid}.XXXXXX)
        picard SamToFastq\
          I={input[0]} F=$tdir/R1.fq F2=$tdir/R2.fq\
          CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 NON_PF=true TMP_DIR=$tdir

        bowtie2 -p {threads} -X 1200\
          -x {config[indexes][bowtie2]}\
          -1 $tdir/R1.fq -2 $tdir/R2.fq\
          --rg-id {wildcards.uuid} --rg PL:ILLUMINA \
          -k 100 --very-sensitive-local --score-min L,0,1.6\
        2> {log[0]} | samtools view -b > {output[0]}

        chmod 600 {output[0]}
        rm -rf $tdir
        '''
