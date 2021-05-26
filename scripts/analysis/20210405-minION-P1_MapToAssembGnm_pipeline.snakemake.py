###################################
#######1k minION
#######
###################################



REP_INDEX = ["BJ13-0001","BJ13-0004","BJ13-0016","BJ13-0017","BJ13-0018","BJ13-0019"]



rule all:
	input:
	##fastq stat
		#expand("/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_Stat/{samp}.fastq.stat.txt",samp=REP_INDEX),
	#bwa map
		#expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/bwa/{samp}/{samp}.sorted.bam",samp=REP_INDEX),
	## map to ref  genome
		expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam",samp=REP_INDEX),
	##bam_index
		expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam.bai",samp=REP_INDEX),
	##bam_depth
		#expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/nglmr_Depth/{samp}_nglmrMapTo_{samp}.sorted.bam.depth.txt",samp=REP_INDEX),	
	##mpileup
		#expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/Mpileup/{samp}_nglmrMapTo_{samp}.sorted.bam.mpileup",samp=REP_INDEX),
	##ivarTrim_varscan2
		#expand("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/varscan2/{samp}.varscan2.snp.vcf",samp=REP_INDEX),
	#

rule TrimReads:
	input:
		"/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_pass/{samp}.fastq",
	output:	
		"/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_pass/{samp}.Trimed.fastq",
	shell:
		"python /shared/liuhj/cao-ncov/script/map_pipeline/ni_trim_Filtfastq.py  -i  {input} -o {output} -m 0  -M  100000  -L 30  -R 30 "   #



rule fastqStat:
	input:
		"/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_pass/fastq_files",
	output:
		"/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_Stat/{samp}.fastq.stat.txt",
	params:
		scriptP="/shared/liuhj/tonglv/ICU_tonglv_scripts"
	shell:
		"python {params.scriptP}/fastq_readsLen_dictrib.readsNum_baseNum.py -i {input} -o {output}"


rule NGLMR_map:
	input:
		"/shared/liuhj/tonglv/process/spades_assembly/assembly_fasta/P1-ONT-Samples/{samp}.fasta",
		"/shared/liuhj/tonglv/data/20210103-ONT-P1/fastq_pass/{samp}.fastq",
	output:	
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam",
	threads:30
	shell:
		"ngmlr -t 30 -r {input[0]}  -q {input[1]} -x ont | samtools sort -@  30 |  samtools view -F 4 -o   {output}"   #-t {nglmr_THREADS}




#
#rule bam_sort:
	#input:
		#"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmr2.sam",
	#output:	
		#temp("/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmr2.sorted.sam"),
	#threads:20
	#shell:
		#"samtools sort   {input}  -o  {output}  "   #-t {nglmr_THREADS}
#
#rule sam_2_bam:
	#input:
		#"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmr2.sorted.sam",
	#output:	
		#"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam",
	#threads:20
	#shell:
		#"samtools view   -o  {output}  {input} "   #-t {nglmr_THREADS}

#
#rule NGLMR_bamsort:
	#input:
		#"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmr2.bam",
	#output:
		#"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam"
	#shell:
		#"samtools sort -@ 6 {input} -o {output}"


rule bam_index:
	input:
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam"
	output:
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam.bai"
	shell:
		"samtools index  {input}"


rule bam_depth:
	input:
		"/shared/liuhj/tonglv/ref_tonglv_AE004091.2/AE004091.2.fna/ref_genome_NC04512.2.fasta",
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam",
	output:
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/nglmr_Depth/{samp}_nglmrMapTo_{samp}.sorted.bam.depth.txt"
	shell:
		"samtools depth -a  -d 30000  --reference {input[0]} {input[1]}  >{output}"


rule mpileup:
	input:
		"/shared/liuhj/tonglv/ref_tonglv_AE004091.2/AE004091.2.fna/ref_genome_NC04512.2.fasta",
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam",
		"/shared/liuhj/tonglv/ref_tonglv_AE004091.2/AE004091.2.fna/ref_genome_NC04512.2.fasta.fai",
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/{samp}_nglmrMapTo_{samp}.sorted.bam.bai"
	output:
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/Mpileup/{samp}_nglmrMapTo_{samp}.sorted.bam.mpileup"
	log:
		"/shared/liuhj/tonglv/process/20210405-ONT-P1-mapTo-SampAssembGnm/Mpileup/{samp}.ivar_trim.mpileup.log"
	shell:
		"samtools mpileup -A  -d 10000  -B -Q 0 --reference  {input[0]}  {input[1]}  1>{output} 2>{log}"     

#-Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [13]
#-q, --min-MQ INT        skip alignments with mapQ smaller than INT [0]
#-d, –max-depth 最大测序深度，过滤掉超深度测序的位点





