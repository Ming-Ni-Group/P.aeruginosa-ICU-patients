###################################
#######
#######
###################################

samps_lst = ["BJ13-0001_BDMS190613907-1a","BJ13-0004_BDMS190602537-1a","BJ13-0005_BDMS190613910-2a","BJ13-0006_BDMS190602539-1a","BJ13-0007_BDMS190624850-1a","BJ13-0008_BDMS190602541-1a","BJ13-0009_BDMS190624851-1a","BJ13-0010_BDMS200010509-1a","BJ13-0011_BDMS190602544-1a","BJ13-0014_BDMS190602547-1a","BJ13-0015_BDMS200010510-1a","BJ13-0016_BDMS190602549-1a","BJ13-0017_BDMS200010530-1a","BJ13-0018_BDMS190602551-1a","BJ13-0019_BDMS200010532-1a","BJ13-0028_BDMS190624856-1a","BJ13-0042_BDMS190602572-1a","BJ13-0044_BDMS200010520-1a","BJ13-0046_BDMS200010521-1a","BJ13-0047_BDMS190602575-1a","BJ13-0048_BDMS200010522-1a","BJ13-0052_BDMS190602580-1a","BJ13-0058_BDMS200010524-1a","BJ13-0059_BDMS200010527-1a","BJ13-0060_BDMS190624867-1a","BJ13-0061_BDMS190602586-1a","BJ13-0064_BDMS200010514-1a","BJ13-0065_BDMS200010528-1a","BJ13-0067_BDMS190602592-1a","BJ13-0068_BDMS190624884-1a","BJ13-0069_BDMS190624889-1a","BJ13-0071_BDMS190624878-1a","BJ13-0072_BDMS190602597-1a","BJ13-0073_BDMS190602598-1a","BJ13-0074_BDMS200010526-1a","BJ13-0075_BDMS200010525-1a","BJ13-0076_BDMS200010513-1a","BJ13-0077_BDMS190602602-1a","BJ13-0078_BDMS200008814-1a","BJ13-0079_BDMS190602604-1a","BJ13-0081_BDMS190624882-1a","BJ13-0082_BDMS190624897-1a","BJ13-0084_BDMS190602609-1a","BJ13-0085_BDMS190602610-1a","BJ13-0086_BDMS200008815-1a","BJ13-0087_BDMS190602612-1a","BJ13-0088_BDMS190624888-1a","BJ13-0089_BDMS190624903-1a","BJ13-0091_BDMS190624891-1a","BJ13-0093_BDMS190613960-1a","BJ13-0097_BDMS190602622-1a","BJ13-0098_BDMS190602623-1a","BJ13-0100_BDMS200010512-1a","BJ13-0101_BDMS190624869-1a","BJ13-0102_BDMS200010529-1a","BJ13-0105_BDMS190602630-1a","BJ13-0108_BDMS190602633-1a","BJ13-0109_BDMS200010533-1a","BJ13-0110_BDMS190602635-1a","BJ13-0114_BDMS190602639-1a","BJ13-0120_BDMS190602645-1a","BJ13-0121_BDMS190624873-1a","BJ13-0122_BDMS200010531-1a","BJ13-0123_BDMS200010511-1a","BJ13-0124_BDMS190602649-1a","BJ13-0127_BDMS190602652-1a","BJ13-0128_BDMS190624904-1a","BJ13-0134_BDMS190602659-1a","BJ13-0140_BDMS200010518-1a","BJ13-0147_BDMS200010519-1a","BJ13-0148_BDMS190624877-1a"]



rule all:
	input:
	## bwa to ref  genome
		expand("/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna.bwt"),
		expand("/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",samp=samps_lst),
		expand("/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam.bai",samp=samps_lst),
	##cnvkit 
		expand("/shared/liuhj/tonglv/process/20210406-cnvkit-call-CNV/{samp}",samp=samps_lst),
		expand("/shared/liuhj/tonglv/process/20210406-cnvkit-call-CNV/reliableGainedCNV/{samp}.sort.rmdup.bintest.reliableGainedCNV.cns",samp=samps_lst)

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
		temp("/shared/liuhj/tonglv/process/bam/{samp}_bwa_{OtherBacteriaRef}.bam")
	shell:
		"bwa  mem  -t  5  {input[0]}  {input[1]}  {input[2]} \
      		|  samtools view -bS -F 4 -  > {output}   "   #

##add picard mark dump



rule bam_index:
	input:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam"
	output:
		"/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam.bai"
	shell:
		"samtools  sort  --threads   5  {input}  -o {output}"



rule cnvkit_call_cnv:
	input:
		NGS_bwa_bam="/shared/liuhj/tonglv/process/bam/{samp}.sort.rmdup.bam",
		Simulate_PAreadsMap_bam="/shared/liuhj/tonglv/process/20210405-PAsimulateReads-mapToAE0049/PA.refCutToEcoli.step1.Len150_bwa_AE004091.2.sort.bam",
		refGnmF="/shared/liuhj/tonglv/ref/ref_tonglv_AE004091.2/AE004091.2.fna"
	output:
		outP=directory("/shared/liuhj/tonglv/process/20210406-cnvkit-call-CNV/{samp}")
	params:
		threads="20",
		target_avg_size="150"
	shell:
		"cnvkit.py  batch   {input.NGS_bwa_bam}  -n  {input.Simulate_PAreadsMap_bam}   -m wgs  --output-dir  {output.outP}  -f {input.refGnmF}  \
		 -p {params.threads}  --drop-low-coverage  --target-avg-size  {params.target_avg_size}"


rule find_CNV_region:
	input:
		cnvkit_Found_CNVF="/shared/liuhj/tonglv/process/20210406-cnvkit-call-CNV/{samp}/{samp}.sort.rmdup.bintest.cns"
	output:
		Filt_reliableCNVRegion="/shared/liuhj/tonglv/process/20210406-cnvkit-call-CNV/reliableGainedCNV/{samp}.sort.rmdup.bintest.reliableGainedCNV.cns"

	shell:
		"cat {input.cnvkit_Found_CNVF} | awk '$6>0 '   > {output.Filt_reliableCNVRegion}"





