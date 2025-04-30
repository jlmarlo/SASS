

def find_files(wildcards):
    ck_output = checkpoints.MakeVariantLists.get(**wildcards).output[0]
    seeds, = glob_wildcards(os.path.join(ck_output, "seed{seed}variantlist.tsv"))
    return expand(os.path.join(ck_output, "{seed}variantlist.tsv"), seed=seeds)



rule MakeTSVs:
    input:
        SNPVCF = config['snpvcf'], 
        IndelVCF = config['indelvcf']
    output:
        SNPTSV = "outputs/TSVFiles/SNPs.tsv",
        IndelTSV = "outputs/TSVFiles/Indels.tsv"
    shell:
        '''
            gatk VariantsToTable \
                -V {input.SNPVCF} \
                -F CHROM -F POS -F REF -F ALT \
                -O {output.SNPTSV}
            gatk VariantsToTable \
                -V {input.IndelVCF} \
                -F CHROM -F POS -F REF -F ALT \
                -O {output.IndelTSV}
        '''
checkpoint MakeVariantLists:
    input:
        SNP = "outputs/TSVFiles/SNPs.tsv",
        Indel = "outputs/TSVFiles/Indels.tsv"
    output:
        directory('outputs/VariantSets')
    params:
        numList = config['simGen'],
        numSnps = config['snps'],
        numIndels = config['indels'],
        seed = config['seed'],
        rLibrary = config['library'],
        prefix = config['output']
    resources:
    shell:
        r'''
            scripts/SelectingVariantsforSimulation.R {params.rLibrary} {input.SNP} {input.Indel} {params.numSnps} {params.numIndels} {output}/{params.prefix} {params.numList} {params.seed}
        '''

rule grabSNPs:
    input:
        SNPList = "outputs/VariantLists/seed{Seed}variantregionlistSNPs.tsv"
    output:
        SNPVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}SNPs.vcf.gz"
    params:
        SnpPoolVCF = config['snpvcf']
    resources:
        time = 360,
        mem_mb = 10000
    shell:
        '''
            bcftools view {params.SnpPoolVCF} \
                -T {input.SNPList} \
                -Oz \
                -o {output.SNPVCF}
            gatk IndexFeatureFile -I {output.SNPVCF}
        '''

rule combineMASNPs:
    input:
        SNPVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}SNPs.vcf.gz"
    output:
        NormVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}.norm.vcf.gz"
    resources:
        time = 360,
        mem_mb = 10000
    shell:
        '''
            bcftools norm {input.SNPVCF} \
                -d snps \
                -m + \
                -Oz \
                -o {output.NormVCF}
            gatk IndexFeatureFile -I {output.NormVCF}
        '''

rule grabIndels:
    input:
        SeedList = "outputs/SeedLists/seed{Seed}variantregionlistindels.tsv"
    output:
        IndelVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}Indels.vcf.gz"
    params:
        IndelsPool = config['indelvcf'] 
    resources:
        time = 180,
        mem_mb = 10000
    shell:
        '''
            bcftools view {params.IndelsPool} \
                -T {input.SeedList} \
                -Oz \
                -o {output.IndelVCF}
            gatk IndexFeatureFile -I {output.IndelVCF}
        '''

rule SitesOnlySNPs:
    input:
        SNPVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}.norm.vcf.gz"    
    output:
        SOVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}SNPs_SO.vcf.gz"
    resources:
        time = 360,
        mem_mb = 10000
    shell:
        '''
            gatk MakeSitesOnlyVcf \
                -I {input.SNPVCF}\
                -O {output.SOVCF}
        '''

rule SitesOnlyIndels:
    input:
        IndelVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}Indels.vcf.gz"
    output:
        SOVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}Indels_SO.vcf.gz"
    resources:
        time = 180,
        mem_mb = 10000
    shell:
        '''
            gatk MakeSitesOnlyVcf \
                -I {input.IndelVCF} \
                -O {output.SOVCF}
        '''

rule CombineVCFs:
    input:
        SnpVCF = "outputs/SeedVCFsIntermediate/Seed{Seed}SNPs_SO.vcf.gz",
        IndelVCF ="outputs/SeedVCFsIntermediate/Seed{Seed}Indels_SO.vcf.gz"
    output:
        Gold = "outputs/SeedVCFs/Seed{Seed}Golden.vcf"
    resources:
        time = 720,
        mem_mb = 10000
    shell:
        '''
            bcftools merge \
                -m both \
                --force-samples \
                {input.SnpVCF} \
                {input.IndelVCF} \
                -Ov -o {output.Gold}
            gatk IndexFeatureFile -I {output.Gold}
        '''

rule BCFStats:
    input:
        Gold = "outputs/SeedVCFs/Seed{Seed}Golden.vcf"
    output:
        Stats = "outputs/BCFStats/Seed{Seed}InsertionStats.txt"
    resources:
        time = 60,
        mem_mb = 5000
    shell:
        '''
            bcftools stats {input.Gold} > {output.Stats}
        '''


def find_files(wildcards):
    ck_output = checkpoints.MakeVariantLists.get(**wildcards).output[0]
    dir = 'outputs/BCFStats/'
    seeds, = glob_wildcards(os.path.join(ck_output, "seed{seed}variantlist.tsv"))
    return expand(os.path.join(dir, "Seed{seed}InsertionStats.txt"), seed=seeds)

rule gather_stats:
    input:
        find_files
    output:
        touch('outputs/BCFStats/all.done')
    shell:
        "echo All done"
