import glob
import os
import re

def find_bams(wildcards):
    sample = wildcards.sample
    print(repr(sample))
    file = sample+'*.*am'
    index = sample+'*.*ai'
    print(file, index)
    dir = "datafiles/bams"
    filePath = []
    filePath = glob.glob(os.path.join(dir,file))
    indexPath = glob.glob(os.path.join(dir,index))
    if not filePath or not indexPath:
        raise ValueError("No bam found for sample:"+sample)
    elif len(filePath) >1 or len(indexPath) > 1:
        raise ValueError("Too many bams found for sample:"+sample)
    bam = filePath[0]
    index = indexPath[0]
    return bam,index

def find_fastqs(wildcards):
    sample = wildcards.sample

    #directory
    dir = "datafiles/fastqs"
    
    ##So far as I know this is all the formats I've seen
    options = [sample+"*R1*.f*q.gz",sample+"*read1*.f*q.gz",sample+"*_1.f*q.gz"]
    r1_files = []
    for option in options:
        r1_files = glob.glob(os.path.join(dir,option))
        if (len(r1_files)) > 0:
            break
    if len(r1_files) ==0:
        raise ValueError(f"No R1 FASTQ found for sample {sample}")
    r1=r1_files[0]
    if 'R1' in r1:
        r2 = r1.replace("R1","R2")
    elif 'read1' in r1:
        r2 = r1.replace("read1","read2")
    elif '_1.f' in r1:
        r2 = r1.replace('_1.f','_2.f')
    else:
        raise ValueError(f"Cannot find R2 file for file: {r1}")
    return r1,r2


rule genomecov:
    input:
        find_bams
    output:
        genomecov = "outputs/Models/genomcov/{sample}.genomecov"
    params:
        reference = config["ref"]
    resources:
        time = 1440,
        mem_mb = 100000
    conda:
        "Python2"
    shell:
        '''
            bedtools genomecov \
                -d \
                -ibam {input[0]} \
                -g {params.reference} > {output.genomecov}
        '''
rule computeGC:
    input:
        genomecov = "outputs/Models/genomcov/{sample}.genomecov"
    output:
        gcmodel = "outputs/Models/GCModels/{sample}GCModel"
    params:
        reference = config["ref"],
        program = 'scripts/neat-genreads/utilities/'
    resources:
        time = 1440,
        mem_mb = 100000
    conda:
        "Python2"
    shell:
        '''
            python {params.program}computeGC.py \
                -r {params.reference} \
                -i {input.genomecov} \
                -w 50 \
                -o {output.gcmodel}
        '''

rule computeFragLength:
    input:
        find_bams
    output: 
        fragmentLength = 'outputs/Models/FragLengthModels/{sample}/{sample}.done',
    resources:
        time = 1440,
        mem_mb = 100000
    params:
        out_dir = 'outputs/Models/FragLengthModels/{wildcards.sample}',
        program = '../../../../scripts/neat_genreads/utilities/'
    conda:
        "Python2"
    shell:
        '''
            mkdir {params.out_dir}
            cd {params.out_dir}
            samtools view {input[0]} | python {params.program}computeFraglen.py
            touch {output.fragmentLength}
        '''

rule computeError:
    input:
        find_fastqs
    output:
        ErrorModel = "outputs/Models/ErrorModels/{sample}ErrorModel"
    params:
        program= 'scripts/neat_genreads/utilities/'
    resources:
        time = 2880,
        mem_mb = 100000
    conda:
        "Python2"
    shell:
        '''
            python {params.program}genSeqErrorModel.py \
                -i {input[0]} \
                -i2 {input[1]} \
                -o {output.ErrorModel}
        '''

