###################################
#######
#######
###################################


REP_INDEX = ["P1-S10","P1-S11","P1-S12","P1-S13","P1-S14","P1-S15","P1-S1","P1-S2","P1-S3","P1-S4","P1-S5","P1-S6","P1-S7","P1-S8","P1-S9","P2-S1","P3-S10","P3-S11","P3-S12","P3-S13","P3-S14","P3-S15","P3-S16","P3-S1","P3-S2","P3-S3","P3-S4","P3-S5","P3-S6","P3-S7","P3-S8","P3-S9","P4-S10","P4-S1","P4-S2","P4-S3","P4-S4","P4-S5","P4-S6","P4-S7","P4-S8","P4-S9","P5-S10","P5-S11","P5-S12","P5-S13","P5-S14","P5-S15","P5-S16","P5-S17","P5-S1","P5-S2","P5-S3","P5-S4","P5-S5","P5-S6","P5-S7","P5-S8","P5-S9","P6-S1","P6-S2","P6-S3","P6-S4","P6-S5","P6-S6","P7-S1","P7-S2","P7-S3","P7-S4","P7-S5","P7-S6"]



rule all:
	input:
		expand("/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.bwt"),
		expand("/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",samp=REP_INDEX),
		expand("/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam.bai",samp=REP_INDEX),
		expand("/shared/liuhj/tonglv/process/mpileup/{samp}.mpileup2Indel.vcf",samp=REP_INDEX),


rule ref_index:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"
	output:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.bwt"
	shell:
		"bwa  index  {input}"


rule bwa_run:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
		"/shared/liuhj/tonglv/data/BGI/data_analysis/{samp}/{samp}_1.clean.fq.gz",
		"/shared/liuhj/tonglv/data/BGI/data_analysis/{samp}/{samp}_2.clean.fq.gz",
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.bwt"
	output:
		temp("/shared/liuhj/tonglv/process/bam/{samp}.bam")
	shell:
		"bwa  mem  -t  5  {input[0]}  {input[1]}  {input[2]} \
      		|  samtools view -bS -F 4 -  > {output}   "   #



rule bam_sort:
	input:
		"/shared/liuhj/tonglv/process/bam/{samp}.bam"
	output:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.bam"
	shell:
		"samtools  sort  --threads   5  {input}  -o {output}"


      	
rule bam_rmdup:
	input:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.bam"
	output:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",
		temp("/shared/liuhj/Ecoli-202104/process/bam/{samp}.sort.rmdup.metrics"),
	shell:
		"picard  MarkDuplicates  INPUT={input}   \
OUTPUT={output[0]}  METRICS_FILE={output[1]}   REMOVE_DUPLICATES=true   "



rule bam_index:
	input:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam"
	output:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam.bai"
	shell:
		"samtools index  {input}"



rule refGonm_faidx:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"
	output:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.fai"
	shell:
		"samtools  faidx  {input}"




rule mpileup:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.fai",
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam.bai"
	output:
		temp("/shared/liuhj/tonglv/process/mpileup/{samp}.mpileup")
	shell:
		"samtools mpileup -A  -d 10000  -B -Q 0 --reference  {input[0]}  {input[1]}  1>{output} "     

#-Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
#-q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
#-d, –max-depth 最大测序深度，过滤掉超深度测序的位点





rule varscan2_mpileup2Indel:
	input:
		"/shared/liuhj/tonglv/process/mpileup/{samp}.mpileup"
	output:
		"/shared/liuhj/tonglv/process/mpileup/{samp}.mpileup2Indel.vcf"
	shell:
		"varscan   pileup2indel  {input}  --min-coverage 50  --min-reads2  10    --min-var-freq  0.1   --variants  indel --output-vcf 1  >{output} "     


#
#Min coverage:	8
#Min reads2:	2
#Min var freq:	0.01
#Min avg qual:	15
#P-value thresh:	0.99
#USAGE: java -jar VarScan.jar pileup2cns [pileup file] OPTIONS
	#pileup file - The SAMtools pileup file
#
	#OPTIONS:
	#--min-coverage	Minimum read depth at a position to make a call [8]
	#--min-reads2	Minimum supporting reads at a position to call variants [2]
	#--min-avg-qual	Minimum base quality at a position to count a read [15]
	#--min-var-freq	Minimum variant allele frequency threshold [0.01]
	#--min-freq-for-hom	Minimum frequency to call homozygote [0.75]
	#--p-value	Default p-value threshold for calling variants [99e-02]
	#--variants	Report only variant (SNP/indel) positions [0]
#



rule snpeffAnno:
	input:
		inputvcfP="/shared/liuhj/tonglv/process/20210323-pindel/deletion-vcfs/"
	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf_snpeff")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		dbID="AE004091",
		snpeffSoftP="/shared/liuhj/software/snpEff",
	shell:
		"bash {params.scriptPath}/pipeline-snpeff.sh  {input} {output} {params.dbID} {params.snpeffSoftP}"

'''
vcfP=$1
outSnpeffVcfP=$2
dbID=$3
snpeffSoftP=$4
'''






