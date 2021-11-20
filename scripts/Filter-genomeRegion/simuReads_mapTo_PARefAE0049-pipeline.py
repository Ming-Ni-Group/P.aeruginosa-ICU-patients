###################################
#######
#######
###################################

REP_INDEX = ["PA.refCutToEcoli.step1.Len150",]
RefGnms = ["AE004091.2"]



rule all:
	input:
		#expand("/shared/liuhj/tonglv/data/BGI/data_analysis/{samp}/{samp}_1.clean.fq",samp=REP_INDEX),
		#expand("/shared/liuhj/tonglv/data/BGI/data_analysis/{samp}/{samp}_2.clean.fq",samp=REP_INDEX),
	## bwa to ref  genome
		expand("/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna.bwt",OtherBacteriaRef=RefGnms),
		expand("/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam",samp=REP_INDEX,OtherBacteriaRef=RefGnms),
		expand("/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam.bai",samp=REP_INDEX,OtherBacteriaRef=RefGnms),
	#filt homo Region
		#expand("/shared/liuhj/tonglv/process/20210402-bwaToOtherBacteria/PA_homoRegion/{samp}_bwa_{OtherBacteriaRef}.PA_homoRegion.txt",samp=REP_INDEX,OtherBacteriaRef=RefGnms),




rule ref_index:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna"
	output:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna.bwt"
	shell:
		"bwa  index  {input}"


rule bwa_run:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna",
		"/shared/liuhj/tonglv/process/20210323-map2_Ecoli/{samp}.fq",
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna.bwt"
	output:
		temp("/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.bam")
	shell:
		"bwa  mem  -t  5  {input[0]}  {input[1]}  \
      		|  samtools view -bS -F 4 -  > {output}   "   #



rule bam_sort:
	input:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.bam"
	output:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam"
	shell:
		"samtools  sort  --threads   5  {input}  -o {output}"

      	
rule bam_rmdup:
	input:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam"
	output:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.rmdup.bam",
		temp("/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.rmdup.metrics"),
	shell:
		"picard  MarkDuplicates  INPUT={input}   \
OUTPUT={output[0]}  METRICS_FILE={output[1]}   REMOVE_DUPLICATES=true   "




rule bam_index:
	input:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam"
	output:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam.bai"
	shell:
		"samtools index  {input}"


rule refGonm_faidx:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna",
	output:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna.fai",
	shell:
		"samtools  faidx  {input}"



rule bam_depth:
	input:
		"/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/{OtherBacteriaRef}.fna",
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam",
	output:
		"/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/depth/{samp}_bwa_{OtherBacteriaRef}.sort.depth.txt"
	shell:
		"samtools depth -a  --reference {input[0]} {input[1]}  >{output}"





rule bam_filtHomoSeqRegion:
	input:
		bamF="/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/{samp}_bwa_{OtherBacteriaRef}.sort.bam",
	output:
		"/shared/liuhj/tonglv/process/20210402-bwaToOtherBacteria/PA_homoRegion/{samp}_bwa_{OtherBacteriaRef}.PA_homoRegion.txt"
	params:
		refID="{OtherBacteriaRef}",
		scriptP="/shared/liuhj/tonglv/ICU_tonglv_scripts"
	shell:
		"python {params.scriptP}/PA-simulateReads-find-homoRegions.py -i {input[0]} -o  {output}  -r {params.refID} "






