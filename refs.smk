#! /usr/bin/env python
# -*- coding: utf-8 -*-
import re
import gzip

localrules: download_reference, extract_genome, extract_transcriptome

rule download_reference:
    output:
        'refs/downloads/{f}'
    params:
        url = lambda wildcards: config['downloads'][wildcards.f]['url'],
        md5 = lambda wildcards: config['downloads'][wildcards.f]['md5']
    shell:
        'mkdir -p $(dirname {output[0]})'
        ' && '
        'curl'
        ' -L'
        ' {params.url}'
        ' > {output[0]}'
        ' && '
        'echo {params.md5}  {output[0]} | md5sum -c -'


rule extract_genome:
    input:
        'refs/downloads/GRCh38.d1.vd1.fa.tar.gz'
    output:
        config['sequences']['genome'],
        config['sequences']['genome_idx'],        
        config['sequences']['genome_dict']
    conda:
        "envs/refutils.yaml"
    shell:
        'mkdir -p $(dirname {output[0]})'
        ' && '
        'tar -Oxzf {input} | bgzip > {output[0]}'
        ' && '
        'samtools faidx {output[0]}'
        ' && '        
        'picard CreateSequenceDictionary R={output[0]} O={output[2]}'


rule extract_transcriptome:
    input:
        'refs/downloads/gencode.v22.annotation.gtf.gz',
        config['sequences']['genome']
    output:
        config['sequences']['transcripts'],
        config['sequences']['transcripts_dupinfo'],        
        config['sequences']['transcripts_list']
    conda:
        "envs/refutils.yaml"
    shell:
        'tfa=$(mktemp -p {config[tmpdir]})'
        ' && '
        'gunzip -c {input[1]} > $tfa'
        ' && '
        'gunzip -c {input[0]} | gffread - -M -d {output[1]} -g $tfa -w {output[0]}'
        ' && '
        'grep ">" {output[0]} | sed "s/^>//" | sort | uniq > {output[2]}'
        ' && '
        'rm -f $tfa*'


rule id_mapping:
    input:
        'refs/downloads/gencode.v22.annotation.gtf.gz',
        config['sequences']['transcripts_list']
    output:
        config['annotations']['ttg'],
        config['annotations']['gsym']
    run:
        tx_g = {}
        g_sym = {}
        
        raw = (bl.decode('utf-8') for bl in gzip.open(input[0], 'rb'))
        lines = (r.strip('\n').split('\t') for r in raw if not r.startswith('#'))
        for l in lines:
            d = dict(t for t in re.findall('(\S+)\s+"([\s\S]*?)";', l[8]))
            if l[2] == 'gene':
                if d['gene_id'] in g_sym:
                    assert g_sym[d['gene_id']] == d['gene_name'], "Gene name mismatch: %s %s" % (d['gene_name'], g_sym[d['gene_id']])
                g_sym[d['gene_id']] = d['gene_name']
            
            if l[2] == 'transcript':
                if d['transcript_id'] in tx_g:
                    assert tx_g[d['transcript_id']] == d['gene_id'], "Gene ID mismatch: %s %s" % (d['gene_id'], tx_g[d['transcript_id']])
                tx_g[d['transcript_id']] = d['gene_id']
        
        # Generate ttg
        with open(output[0], 'w') as outh:
            print('TXNAME\tGENEID', file=outh)
            txlist = (l.strip('\n').split()[0] for l in open(input[1], 'rU'))
            for tx in txlist:
                print('%s\t%s' % (tx, tx_g[tx]), file=outh)
        
        # Generate gsym
        with open(output[1], 'w') as outh:
            print('GENEID\tSYM', file=outh)
            for t in g_sym.items():
                print('%s\t%s' % t, file=outh)


localrules: telescope_annotation, make_herv_tsv, make_l1_tsv

rule telescope_annotation:
    input:
        config['annotations']['herv'],
        config['annotations']['l1'],
    output:
        config['annotations']['retro']
    shell:
        'mkdir -p $(dirname {output[0]})'
        ' && '
        'cat {input[0]} {input[1]}'
        ' | grep -v "^#"'
        ' | perl -lane \'print if $F[2]=~/exon/\''
        ' | python scripts/sortgtf.py'
        ' > {output[0]}'


rule make_herv_tsv:
    input:
        config['annotations']['herv']
    output:
        config['annotations']['herv_tsv']
    run:
        fields = ['locus', 'class', 'family', 'category', 'chrom', 'start', 'end', 'strand']
        lines = (l.strip('\n').split('\t') for l in open(input[0], 'r') if not l.startswith('#'))
        with open(output[0], 'w') as outh:
            print('\t'.join(fields), file=outh)
            for l in lines:
                if l[2] == 'gene':
                    d = dict(t for t in re.findall('(\S+)\s+"([\s\S]*?)";', l[8]))
                    d['chrom'] = l[0]
                    d['start'] = l[3]
                    d['end'] = l[4]
                    d['strand'] = l[6]
                    d['class'] = "HERV"
                    # d['family'] = d['intModel']
                    d['family'] = d['locus'].split('_')[0]
                    print('\t'.join(d[f] for f in fields), file=outh)


rule make_l1_tsv:
    input:
        config['annotations']['l1']
    output:
        config['annotations']['l1_tsv']
    run:
        fields = ['locus', 'class', 'family', 'category', 'chrom', 'start', 'end', 'strand']
        lines = (l.strip('\n').split('\t') for l in open(input[0], 'r') if not l.startswith('#'))
        with open(output[0], 'w') as outh:
            print('\t'.join(fields), file=outh)
            for l in lines:
                if l[2] == 'exon':
                    d = dict(t for t in re.findall('(\S+)\s+"([\s\S]*?)";', l[8]))
                    d['chrom'] = l[0]
                    d['start'] = l[3]
                    d['end'] = l[4]
                    d['strand'] = l[6]
                    d['class'] = "L1"
                    d['family'] = "L1"
                    print('\t'.join(d[f] for f in fields), file=outh)


localrules: kallisto_index, bowtie2_index

rule kallisto_index:
    input:
        config['sequences']['transcripts']
    output:
        config['indexes']['kallisto']
    conda:
        "envs/kallisto.yaml"
    shell:
        'mkdir -p $(dirname {output})'
        ' && '
        'kallisto index -i {output} {input}'


rule bowtie2_index:
    input:
        config['sequences']['genome']
    output:
        expand(config['indexes']['bowtie2'] + ".{i}.bt2", i = range(1, 5)),
        expand(config['indexes']['bowtie2'] + ".rev.{i}.bt2", i = range(1, 3))        
    conda:
        "envs/bowtie2.yaml"
    threads: snakemake.utils.available_cpu_count() 
    shell:
        'mkdir -p $(dirname {config[indexes][bowtie2]})'
        ' && '
        'tfa=$(mktemp -p {config[tmpdir]})'
        ' && '
        'gunzip -c {input[0]} > $tfa'
        ' && '
        'bowtie2-build --threads {threads} $tfa {config[indexes][bowtie2]}'
        ' && '
        'rm -f $tfa'


localrules: gsea_download

rule gsea_download:
    output:
        'refs/annotation/{f}.gmt'
    shell:
        'mkdir -p refs/annotation'
        ' && '
        'wget -O {output} https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.2/{wildcards.f}.gmt'

