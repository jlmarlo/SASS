rule genomecov:
	input:
		bam = "/scratch.global/marlo072/Simulating90/bams/{horse}.goldenPath.bam",
		bai = "/scratch.global/marlo072/Simulating90/bams/{horse}.goldenPath.bai"
	output:
		genomecov = "/scratch.global/marlo072/Simulating90/Models/genomcov/{horse}.genomecov"
	params:
		reference = "~/Documents/ReferenceFiles/goldenPath.Ec_build-3.0_wMSY.fa"
	resources:
		time = 1440,
		mem_mb = 100000
	conda:
		"Python2"
	shell:
		'''
			cd ~/neat-genreads/utilities
			bedtools genomecov \
				-d \
				-ibam {input.bam} \
				-g {params.reference}\
				> {output.genomecov}
		'''
rule computeGC:
	input:
		genomecov = "/scratch.global/marlo072/Simulating90/Models/genomcov/{horse}.genomecov"
	output:
		gcmodel = "/scratch.global/marlo072/Simulating90/Models/GCModels/{horse}GCModel"
	params:
		reference = '~/Documents/ReferenceFiles/goldenPath.Ec_build-3.0_wMSY.fa'
	resources:
		time = 1440,
		mem_mb = 100000
	conda:
		"Python2"
	shell:
		'''
			cd ~/neat-genreads/utilities
			python computeGC.py \
				-r {params.reference} \
				-i {input.genomecov} \
				-w 50 \
				-o {output.gcmodel}
		'''
#mkdir /scratch.global/marlo072/Simulating90/Models/FragLengthModels/{wildcards.horse}
#touch {output.start}
rule computeFragLength:
	input:
		bam = "/scratch.global/marlo072/Simulating90/bams/{horse}.goldenPath.bam",
		bai = "/scratch.global/marlo072/Simulating90/bams/{horse}.goldenPath.bai"
	output: 
		fragmentLength = '/scratch.global/marlo072/Simulating90/Models/FragLengthModels/{horse}/{horse}.done',
		#start = '/scratch.global/marlo072/Simulating90/Models/FragLengthModels/{horse}.start'
	resources:
		time = 1440,
		mem_mb = 100000
	conda:
		"Python2"
	shell:
		'''
			cd /scratch.global/marlo072/Simulating90/Models/FragLengthModels/{wildcards.horse}
			samtools view {input.bam} | python ~/neat-genreads/utilities/computeFraglen.py
			touch {output.fragmentLength}
		'''

rule computeError:
	input:
		fastq1 = "/scratch.global/marlo072/Simulating90/fastqs/{horse}_R1.fastq.gz",
		fastq2 = "/scratch.global/marlo072/Simulating90/fastqs/{horse}_R2.fastq.gz"
	output:
		ErrorModel = "/scratch.global/marlo072/Simulating90/Models/ErrorModels/{horse}ErrorModel"
	resources:
		time = 2880,
		mem_mb = 100000
	conda:
		"Python2"
	shell:
		'''
			cd ~/neat-genreads/utilities
			python genSeqErrorModel.py \
				-i {input.fastq1} \
				-i2 {input.fastq2} \
				-o {output.ErrorModel}
		'''

