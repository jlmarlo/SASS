rule grabSNPs:
	input:
		SNPList = "/scratch.global/marlo072/Simulating90/SeedLists/seed{Seed}variantregionlistSNPs.tsv"
	output:
		SNPVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}SNPs.vcf.gz"
	params:
		RobVCF = "~/Documents/ReferenceFiles/EquCab3_FILTERED.Equeess_cab_nucl_wChrUn.vcf.gz"
	resources:
		time = 360,
		mem_mb = 10000
	shell:
		'''
			bcftools view {params.RobVCF} \
				-T {input.SNPList} \
				-Oz \
				-o {output.SNPVCF}
			gatk IndexFeatureFile -I {output.SNPVCF}
		'''

rule combineMASNPs:
	input:
		SNPVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}SNPs.vcf.gz"
	output:
		NormVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}.norm.vcf.gz"
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
		SeedList = "/scratch.global/marlo072/Simulating90/SeedLists/seed{Seed}variantregionlistindels.tsv"
	output:
		IndelVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}Indels.vcf.gz"
	params:
		IndelsPool = "~/Documents/ReferenceFiles/prevalent1-60Indels_joint_genotype_indelsfr_20230726.goldenPath.vcf.gz" 
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
		SNPVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}.norm.vcf.gz"	
	output:
		SOVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}SNPs_SO.vcf.gz"
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
		IndelVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}Indels.vcf.gz"
	output:
		SOVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}Indels_SO.vcf.gz"
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
		SnpVCF = "/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}SNPs_SO.vcf.gz",
		IndelVCF ="/scratch.global/marlo072/Simulating90/SeedVCFsIntermediate/Seed{Seed}Indels_SO.vcf.gz"
	output:
		Gold = "/scratch.global/marlo072/Simulating90/SeedVCFs/Seed{Seed}Golden.vcf"
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
		Gold = "/scratch.global/marlo072/Simulating90/SeedVCFs/Seed{Seed}Golden.vcf"
	output:
		Stats = "/scratch.global/marlo072/Simulating90/BCFStats/Seed{Seed}InsertionStats.txt"
	resources:
		time = 60,
		mem_mb = 5000
	shell:
		'''
			bcftools stats {input.Gold} > {output.Stats}
		'''
