###################################
#######1k minION
#######
###################################

personLst = ["P1","P3","P4","P6","P7","P5_withP2","P5-1","P5-2","P5-3","P3-1","P3-2","P3-3"]

#personLst = ["P1"]


rule all:
	input:
	#bwa map
		## iSNV calling
		expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",person=personLst),
		##vcf
		expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf",person=personLst),
		##snpeff vcf
		expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf_snpeff",person=personLst),

		##general snvFreq_distribut_Stat
		expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/general-snvFreq_distribut_Stat",person=personLst),

		## snvNum_Stat
		expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/snvNum_Stat_noSNV_OtherBacfilted-FiltSimutPARepeat",person= personLst),

		##all sa
		#expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/20210331-snvSNP_AnnoStat_OtherBacFilted_RepeatFilted",person=personLst),

		##SNP genome sub SNP posi
		#expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/subGnm",person=personLst),
		#expand("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/subGnm_prokka_gff",person=personLst),




rule ntfreq_2_FreqBigtable:
	input:
		configF="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}_ntfreq_file_list.txt",  ##format in the config file: path/sampID.ntfreq"\n"
		refF="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
	output:	
		outputF="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
		outputF2="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/out.samps-ntfreq.bigtable.txt",

	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		outputP="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}",
		DEP_THRES=50,       ##最低深度
		MIN_DEP_THRES=5,    ##支持alt的最小reads数
		VALID_SIZE_THRES=100,  ##样本有效需要的有效位点数
		FREQ_THRES=0.02,    ###alt的频率最小值
		STRANDED_RATIO_THRES=0.0,     ## reads链偏性检测


	shell:
		"python {params.scriptPath}/step6_read_ntfreq.py  -i  {input.configF} -r {input.refF}  -o {params.outputP}  \
		-D  {params.DEP_THRES} -A {params.MIN_DEP_THRES} -N {params.VALID_SIZE_THRES} -F {params.FREQ_THRES} -S {params.STRANDED_RATIO_THRES} "   #




rule iSNV_freq_distribut:
	input:
		"/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/general-snvFreq_distribut_Stat")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.05,
		maxFreq	= 1.0,
		step = 2,

	shell:
		"python  {params.scriptPath}/iSNVpy_Freq_distribut_Stat.py  -i {input} -o  {output.outputP} -m {params.minFreq}  -M {params.maxFreq}  -s {params.step}  "




rule convert_iSNVpytable_To_vcf :
	input:
		iSNVtable="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
		ntfreqTable = "/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/out.samps-ntfreq.bigtable.txt",
		RefGnmF = "/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"

	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		MinFreq	= 0.05,
		#snv_MaxFreq	= 0.95,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_iSNVtable_To_vcf.py   -i {input.iSNVtable}  \
		-n {input.ntfreqTable}  -r  {input.RefGnmF} -o  {output.outputP} -m {params.MinFreq}  "


##snpeff ann

#for file in ../../process/person_ntfreq_tables/P*/vcf/*.vcf ;do echo $file ;filename=${file##*/};fileID=${filename%.vcf*}; snpEff  ann   -v  AE004091  $file > ${file%/*}/$fileID.snpeff.vcf ;done

rule snpeffAnno:
	input:
		inputvcfP="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf"
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




'''

rule iSNV_snvNumCount:
	input:
		iSNVtable="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
		infoF = "/shared/liuhj/tonglv/Infos/20210311-info.txt",
	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/snvNum_Stat")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.1,
		maxFreq	= 0.9,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_Freq_distribut_Stat_clade.py  -i {input.iSNVtable}  \
		-f {input.infoF}  -o  {output.outputP} -m {params.minFreq}  -M {params.maxFreq}  "

'''




rule iSNV_snvNumCount_noSNV:
	input:
		iSNVtable="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
		infoF = "/shared/liuhj/tonglv/Infos/20210311-info.txt",
		filtHomoPosi="/shared/liuhj/tonglv/process/20210402-bwaToOtherBacteria/PA-HomoRegions-FromOtherBacteria.txt",
		snpEffPath="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/vcf_snpeff",
		RepeatMaskedGnm="/shared/liuhj/tonglv/process/RefGnmRepeatRegion/repeatmasker_out/AE004091.2.fasta.masked.fasta",
		RepeatPosiF="/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/PA.RepeatRegion-Cov30-identy60.Posi.txt"
	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/snvNum_Stat_noSNV_OtherBacfilted-FiltSimutPARepeat")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		#iSNV_FreqThres = 0.05,     ## for examp: iSNV_FreqThres = 0.05, mean iSNV freq :>=0.05 && <0.95  ;SNP: >=0.95
		minFreq	= 0.1,
		maxFreq	= 0.9,
		
	shell:
		"python  {params.scriptPath}/iSNVpy_snvSNP_Stat_clade.py   -i {input.iSNVtable}  \
		-f {input.infoF} -e {input.snpEffPath} -R {input.RepeatPosiF} -o  {output.outputP}  -m {params.minFreq}  -M {params.maxFreq}  -E {input.filtHomoPosi}   "



rule iSNV_SNP_consensusGnm:
	input:
		iSNVtable="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/all.iSNV_with_SNP.pyResults.txt",
		ntfreqTable = "/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/out.samps-ntfreq.bigtable.txt",
		RefGnmF = "/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"

	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/subGnm")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		MinFreq	= 0.9,
		#snv_MaxFreq	= 0.95,
		

	shell:
		"python  {params.scriptPath}/iSNVpy_iSNVtable_RefGnmtihuanSNP.py   -i {input.iSNVtable}  \
		-n {input.ntfreqTable}  -r  {input.RefGnmF} -o  {output.outputP} -m {params.MinFreq}  "






rule consensusGnm_prokkaAnno:
	input:
		inputP="/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/subGnm"
	output:
		outputP=directory("/shared/liuhj/tonglv/process/20210401-person_ntfreq_tables/{person}/subGnm_prokka_gff")
	params:
		scriptPath="/shared/liuhj/tonglv/ICU_tonglv_scripts/iSNV_calling_bwa2_AE004091.2",
		genbankdb="/home/amax/anaconda3/envs/biotools_lhj/db/Pseudomonas/AE004091.2.gb"

	shell:
		"{params.scriptPath}/pipeline-prokka.sh  {input.inputP}  {output.outputP}  {params.genbankdb} "


#gffP=$1
#outP=$2
#genbankdb=$3  ## tonglv --- /home/amax/anaconda3/envs/biotools_lhj/db/Pseudomonas/AE004091.2.gb











