###################################
###############
#######
###################################

##共71个样本
REP_INDEX = ["BJ13-0001_BDMS190613907-1a","BJ13-0004_BDMS190602537-1a",\
"BJ13-0005_BDMS190613910-2a","BJ13-0006_BDMS190602539-1a","BJ13-0007_BDMS190624850-1a",\
"BJ13-0008_BDMS190602541-1a","BJ13-0009_BDMS190624851-1a","BJ13-0010_BDMS200010509-1a",\
"BJ13-0011_BDMS190602544-1a","BJ13-0014_BDMS190602547-1a","BJ13-0015_BDMS200010510-1a",\
"BJ13-0016_BDMS190602549-1a","BJ13-0017_BDMS200010530-1a","BJ13-0018_BDMS190602551-1a",\
"BJ13-0019_BDMS200010532-1a","BJ13-0028_BDMS190624856-1a","BJ13-0042_BDMS190602572-1a",\
"BJ13-0044_BDMS200010520-1a","BJ13-0046_BDMS200010521-1a","BJ13-0047_BDMS190602575-1a",\
"BJ13-0048_BDMS200010522-1a","BJ13-0052_BDMS190602580-1a",\
"BJ13-0058_BDMS200010524-1a","BJ13-0059_BDMS200010527-1a","BJ13-0060_BDMS190624867-1a",\
"BJ13-0061_BDMS190602586-1a","BJ13-0064_BDMS200010514-1a","BJ13-0065_BDMS200010528-1a",\
"BJ13-0067_BDMS190602592-1a","BJ13-0068_BDMS190624884-1a","BJ13-0069_BDMS190624889-1a",\
"BJ13-0071_BDMS190624878-1a","BJ13-0072_BDMS190602597-1a","BJ13-0073_BDMS190602598-1a",\
"BJ13-0074_BDMS200010526-1a","BJ13-0075_BDMS200010525-1a","BJ13-0076_BDMS200010513-1a",\
"BJ13-0077_BDMS190602602-1a","BJ13-0078_BDMS200008814-1a","BJ13-0079_BDMS190602604-1a",\
"BJ13-0081_BDMS190624882-1a","BJ13-0082_BDMS190624897-1a","BJ13-0084_BDMS190602609-1a",\
"BJ13-0085_BDMS190602610-1a","BJ13-0086_BDMS200008815-1a","BJ13-0087_BDMS190602612-1a",\
"BJ13-0088_BDMS190624888-1a","BJ13-0089_BDMS190624903-1a","BJ13-0091_BDMS190624891-1a",\
"BJ13-0093_BDMS190613960-1a","BJ13-0097_BDMS190602622-1a","BJ13-0098_BDMS190602623-1a",\
"BJ13-0100_BDMS200010512-1a","BJ13-0101_BDMS190624869-1a","BJ13-0102_BDMS200010529-1a",\
"BJ13-0105_BDMS190602630-1a","BJ13-0108_BDMS190602633-1a","BJ13-0109_BDMS200010533-1a",\
"BJ13-0110_BDMS190602635-1a","BJ13-0114_BDMS190602639-1a","BJ13-0120_BDMS190602645-1a",\
"BJ13-0121_BDMS190624873-1a","BJ13-0122_BDMS200010531-1a","BJ13-0123_BDMS200010511-1a",\
"BJ13-0124_BDMS190602649-1a","BJ13-0127_BDMS190602652-1a","BJ13-0128_BDMS190624904-1a",\
"BJ13-0134_BDMS190602659-1a","BJ13-0147_BDMS200010519-1a",\
"BJ13-0148_BDMS190624877-1a","BJ13-0140_BDMS200010518-1a",]



#REP_INDEX = ["BJ13-0006_BDMS190602539-1a","BJ13-0007_BDMS190624850-1a"]


rule all:
	input:
	##pindel
		expand("/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_D.vcf",samp=REP_INDEX),  #pindel_indel_F=
	##delly
		expand("/shared/liuhj/tonglv/process/20210323-delly-SV-calling/vcf/{samp}.sort.rmdup.delly-SV.vcf",samp=REP_INDEX)
	




rule pindel_indel:
	input:
		ref_F="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
		bam_F="/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",
	output:	
		configF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.confg.txt",  ##pindel_confg_F=
		delF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_D",  #pindel_indel_F=
		insertF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_SI",  #pindel_indel_F=
	params:
		sampID="{samp}",
		chrom_NAME="AE004091.2",
		pindel_outdir="/shared/liuhj/tonglv/process/20210323-pindel/{samp}",
		min_mapping_Q=30,
		min_alt_reads=10,
		threads=2,
		insertSize=360
		#scriptP="/shared/liuhj/tonglv/ICU_tonglv_scripts"
	log:
		"/shared/liuhj/tonglv/process/20210323-pindel/{samp}.pindel.Indel.log"
	shell:
		#"bash  {params.scriptP}/pindel.sh  {input[0]} {input[1]} {output[0]} {output[1]} {params[0]} {params[1]} {params[2]} {params[3]} {params[4]} {params[5]}  2>{log} "
		"echo -e  {input.bam_F}\"\t\"\"{params.insertSize}\"\"\t\"{params.sampID}   >{output.configF};\
		pindel  -f  {input.ref_F}  -i  {output.configF}  -c ALL -A  {params.min_mapping_Q}  -M  {params.min_alt_reads} -T  {params.threads}  -o  {params.pindel_outdir}/{params.sampID}.pindel"

#-A  the minimal mapping quality of the reads Pindel uses as anchor If you only need high confident calls, set to 30 or higher(default 0)
##-M  --minimum_support_for_event
#-T  threads


rule pindeF_Tovcf:
	input:
		ref_F="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
		DelF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_D",  #pindel_indel_F=
		insertF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_SI",  #pindel_indel_F=

	output:
		defVcfF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_D.vcf",  #pindel_indel_F=
		instVcfF="/shared/liuhj/tonglv/process/20210323-pindel/{samp}/{samp}.pindel_SI.vcf",  #pindel_indel_F=
	params:
		sampID="{samp}",
		chrom_NAME="AE004091.2",
		pindel_outdir="/shared/liuhj/tonglv/process/20210323-pindel/{samp}",
		min_alt_reads=10,

	shell:
		"pindel2vcf -G  -r {input.ref_F}  -R {params.chrom_NAME} -p  {input.DelF}  -d 20210324  -e  {params.min_alt_reads}  -v  {output.defVcfF};\
		pindel2vcf -G  -r {input.ref_F}  -R {params.chrom_NAME} -p  {input.insertF}  -d 20210324  -e  {params.min_alt_reads}  -v  {output.instVcfF};"

##-e  可以设置支持alt的最小reads数



rule SVcall_delly:
	input:
		ref_F="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna",
		bamF="/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",
	output:	
		outbcfF=temp("/shared/liuhj/tonglv/process/20210323-delly-SV-calling/{samp}.sort.rmdup.delly-SV.bcf"),
		outbcf_csiF=temp("/shared/liuhj/tonglv/process/20210323-delly-SV-calling/{samp}.sort.rmdup.delly-SV.bcf.csi"),
	shell:
		"delly  call  -t DEL  -g  {input.ref_F}  -o {output.outbcfF}    {input.bamF}"



rule bcfTovcf:
	input:
		bcfF="/shared/liuhj/tonglv/process/20210323-delly-SV-calling/{samp}.sort.rmdup.delly-SV.bcf",
	output:	
		vcfF="/shared/liuhj/tonglv/process/20210323-delly-SV-calling/vcf/{samp}.sort.rmdup.delly-SV.vcf"
	shell:
		"bcftools convert  -o  {output.vcfF}   {input.bcfF}"








