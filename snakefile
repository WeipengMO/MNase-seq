'''
Author       : windz
Date         : 2020-09-03 13:16:08
LastEditTime : 2020-09-03 22:19:29
Description  : 
'''
import os
from glob import glob


# create log dir
if not os.path.exists('log'):
    os.mkdir('log')

SAMPLE_NAME = [os.path.basename(fn).split('.')[0] for fn in glob('raw_data/*.1.fastq.gz')]

rule all:
    input:
        expand('aligned_data/{sample_name}/{sample_name}.sorted.bam', sample_name=SAMPLE_NAME),
        expand('danpos_result/{sample_name}.o', sample_name=SAMPLE_NAME),
        expand('danpos_result/pooled/aligned_data_{sample_name}.smooth.bw', sample_name=SAMPLE_NAME),


rule run_fastp:
    input:
        'raw_data/{sample_name}.1.fastq.gz',
    output:
        'raw_data/{sample_name}.clean.fastq.gz'
    threads: 8
    shell:
        '''
fastp -w {threads} -i {input} -o {output}
        '''

# 👉比对
rule run_bowtie2:
    input:
        'raw_data/{sample_name}/{sample_name}.clean.fastq.gz'
    output:
        'aligned_data/{sample_name}/{sample_name}.sorted.bam'
    params:
        genome='~/db/Arabidopsis_thaliana/dna/bowtie2_index/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'
    threads: 60
    shell:
        '''
bowtie2 -q -p {threads} -t --no-mixed --no-discordant --no-unal --dovetail -x {params.genome} -U {input} | samtools sort -@ 64 -O bam -o {output} -
samtools index -@ {threads} {output}
        '''


rule run_danpos:
    input:
        dir='aligned_data/{sample_name}/',
        bam='aligned_data/{sample_name}/{sample_name}.sorted.bam'
    output:
        out='danpos_result/{sample_name}.o'
    threads: 1
    shell:
        '''
export PATH="/public/home/mowp/miniconda3/envs/R/bin/:$PATH"
python DANPOS3/danpos.py dpos {input.dir} --out danpos
touch {output.out}
        '''

rule run_danpos_profiles:
    input:
        'danpos_result/pooled/aligned_data_{sample_name}.smooth.wig'
    output:
        'danpos_profile/{sample_name}/{sample_name}.out'
    threads: 1
    params:
        out_dir='danpos_profile/{sample_name}',
        gene_file='supplementary_data/araport11.gene_file.tsv',
        danpos='/public/home/mowp/test/mnase_seq_nucmap/DANPOS3/danpos.py'
    shell:
        '''
export PATH="/public/home/mowp/miniconda3/envs/R/bin/:$PATH"
cd {params.out_dir}
python {params.danpos} profile ../../{input} --genefile_paths ../../{params.gene_file}
touch ../../{output}
        '''

rule wigToBigWig:
    input:
        'danpos_result/pooled/aligned_data_{sample_name}.smooth.wig'
    output:
        'danpos_result/pooled/aligned_data_{sample_name}.smooth.bw'
    params:
        'supplementary_data/Arabidopsis_thaliana.TAIR10.chrom.sizes'
    threads: 1
    shell:
        '''
wigToBigWig {input} {params} {output}
        '''

