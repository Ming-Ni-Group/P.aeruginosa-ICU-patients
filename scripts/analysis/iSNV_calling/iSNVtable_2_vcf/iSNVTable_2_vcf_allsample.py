#!/usr/bin/python
# -*- coding: utf-8 -*-

###### filt kela naiyao weidian reads  and chars  from ONT bam file
import sys
import linecache
import os


###############

def vcf_Header(SAMPID):
	HeaderL = []
	HeaderL.append("##fileformat=VCFv4.1")
	HeaderL.append("##source=iSNV")
	HeaderL.append('##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 20">')
	HeaderL.append('##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">')
	HeaderL.append('##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">')
	HeaderL.append('##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">')
	HeaderL.append('##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">')
	HeaderL.append(str('##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">'))
	HeaderL.append('##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">')
	HeaderL.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
	HeaderL.append('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
	HeaderL.append('##FORMAT=<ID=MDP,Number=1,Type=Integer,Description="Raw Read Depth from mpileup">')
	HeaderL.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= MapQThreds">')
	HeaderL.append('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">')
	HeaderL.append('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">')
	HeaderL.append('##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">')
	#HeaderL.append("##FORMAT=<ID=PVAL,Number=1,Type=String,Description='P-value from Fisher's Exact Test>")
	#HeaderL.append('##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">')
	#HeaderL.append('##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">')
	#HeaderL.append('##FORMAT=<ID=RDF,Number=1,Type=Integer,Description="Depth of reference-supporting bases on forward strand (reads1plus)">')
	#HeaderL.append('##FORMAT=<ID=RDR,Number=1,Type=Integer,Description="Depth of reference-supporting bases on reverse strand (reads1minus)">')
	#HeaderL.append('##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Depth of variant-supporting bases on forward strand (reads2plus)">')
	#HeaderL.append('##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Depth of variant-supporting bases on reverse strand (reads2minus)">')
	HeaderL.append('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	' + SAMPID)

	HEADER = "\n".join(HeaderL)
	return HEADER


def sample_ID(iSNVTable):
	iSNVTablels = open(iSNVTable,'r').readlines()
	##find sample col  index
	for iSNVTablel in iSNVTablels:
		if "pos	iSNV/SNP#" in iSNVTablel:
			sampIDs = iSNVTablel.split("\t")[2:]
			sampIDs.remove("\n")
			
	return sampIDs	



def iSNVsampPos(iSNVTable,SAMPID):
	iSNVTablels = open(iSNVTable,'r').readlines()
	##find sample col  index
	for iSNVTablel in iSNVTablels:
		if "pos	iSNV/SNP#" in iSNVTablel:
			sampIDs = iSNVTablel.split("\t")
			sampIDcount = 0
			for sampID in sampIDs:
				if sampID == SAMPID:
					sampIDindex = sampIDcount
					break
				sampIDcount += 1
			break	

	###sample  posi  that isn't  NO or NA
	sampPosL = []
	sampPosFreqD = {}
 	for iSNVTablel in iSNVTablels:
		if "#" not in iSNVTablel and iSNVTablel != "\n":
		#if "pos	iSNV#" in iSNVTablel:
			sampsFreq = iSNVTablel.split("\t")
			POS = sampsFreq[0]
			sampFREQ = sampsFreq[sampIDindex]
			if sampFREQ != "NA" and sampFREQ != "NO" and float(sampFREQ) <= FREQ_threshold:
				#print POS + "\t" + sampFREQ
				sampPosL.append(POS)
			if  sampFREQ != "NA" and sampFREQ != "NO" and  float(sampFREQ) >= 0.9:	
				sampPosFreqD[POS] = sampFREQ			
	return sampPosL


def samp_posi_info(infoF,sampPosL,SAMPID):
	infoFls = open(infoF,'r').readlines()
#	SNPinfoFls = open(SNPinfoF,'r').readlines()
	sampSiteInfoD = {}
	notHegeL = []
	for infoFl in infoFls:
		if "#" not in infoFl and infoFl != "\n":
			infos = infoFl.split("\t")
			ID = infos[0]
			Site = infos[1]	
			if ID == SAMPID  and Site in sampPosL:	
				REF = infos[4].split(":")[0]
				nt_pattern = infos[18].split(";")
				Anum=int(nt_pattern[0].split(":")[1])
				Gnum=int(nt_pattern[1].split(":")[1])
				Cnum=int(nt_pattern[2].split(":")[1])
				Tnum=int(nt_pattern[3].split(":")[1])
				totnum=int(nt_pattern[4].split(":")[1])
				charsCount = 0
				NumD = {"A":Anum,"G":Gnum,"C":Cnum,"T":Tnum}
				FreqD = {}
				Afreq = str(round(float(Anum)*100/totnum,3))
				FreqD["A"] = Afreq
				Gfreq = str(round(float(Gnum)*100/totnum,3))
				FreqD["G"] = Gfreq
				Cfreq = str(round(float(Cnum)*100/totnum,3))
				FreqD["C"] = Cfreq
				Tfreq = str(round(float(Tnum)*100/totnum,3))
				FreqD["T"] = Tfreq
					
				if Anum > baseCovThreshold:	
					if REF != "A":
						charsCount += 1
				if Gnum > baseCovThreshold:
					if REF != "G":
						charsCount += 1
				if Cnum > baseCovThreshold:
					if REF != "C":
						charsCount += 1
				if Tnum > baseCovThreshold:
					if REF != "T":
						charsCount += 1


				#if Afreq > baseFreqThreds:	
					#if REF != "A":
						#charsCount += 1
				#if Gfreq > baseFreqThreds:
					#if REF != "G":
						#charsCount += 1
				#if Cfreq > baseFreqThreds:
					#if REF != "C":
						#charsCount += 1
				#if Tfreq > baseFreqThreds:
					#if REF != "T":
						#charsCount += 1


				if charsCount <= 1:   ##只选择那些有两种碱基的位点,除了REf，只有一种突变碱基,用baseCovThreshold限制
				#if charsCount >= 1:   ##只选择那些有两种碱基的位点,除了REf，只有一种突变碱基
					MDP = infos[5]   #cov_from_mpileup
					DP = infos[7]   #cov_noIndel_QC
					minor_cov = infos[9]
					second_cov = infos[11]	
					minor_freq = infos[10]
					second_freq = infos[12]						
					chars_max = infos[4].split(":")[1].split("-")[0]
					chars_second = infos[4].split(":")[1].split("-")[1]
					#print REF + "\t" + chars_max + "\t" + chars_second
					if chars_max == REF:
						ALT = chars_second
					else:
						ALT = chars_max

					if chars_second == REF:
						ALT = chars_max

					RD = str(NumD[REF])
					AD = str(NumD[ALT])
					#print RD + "\t"  +AD
					FREQ = FreqD[ALT]

					sampSiteInfoL = []
					sampSiteInfoL.append(REF)
					sampSiteInfoL.append(ALT)
					sampSiteInfoL.append(MDP)
					sampSiteInfoL.append(DP)
					sampSiteInfoL.append(RD)
					sampSiteInfoL.append(AD)
					sampSiteInfoL.append(FREQ)
					sampSiteInfo = "\t".join(sampSiteInfoL)
					#print sampSiteInfo
					sampSiteInfoD[Site] = sampSiteInfo

				else:
					notHegeL.append(Site)
					#print Site + "\t" + REF + "\t" + "Not  hege"
					#print nt_pattern
					#print
	#print "buhegePos:   " + str(len(notHegeL))
	return sampSiteInfoD				



def vcfFomat(sampPosL,snvInfoD,SNPInfoD):
	outL = []
	NotHegeCount = 0
	for sampPos in sampPosL:
		lineL = ''
		if sampPos in snvInfoD.keys():
			lineL = snvInfoD[sampPos].split("\t")
		elif sampPos in SNPInfoD.keys():
			lineL = SNPInfoD[sampPos].split("\t")
		else:
			NotHegeCount += 1

		
		outlineL = []
		if lineL != '':
			outlineL.append(RefChromID)
			outlineL.append(sampPos)
			outlineL.append(".")
			outlineL.append(lineL[0])  ##REF
			outlineL.append(lineL[1])  ##ALT
			outlineL.append(".")
			outlineL.append("PASS")
			outlineL.append(";".join(["ADP="+lineL[3],"WT=0","HET=0","HOM=1","NC=0"]))
			outlineL.append("GT:GQ:MDP:DP:RD:AD:FREQ")
			outlineL.append(":".join(["1/1","255",lineL[2],lineL[3],lineL[4],lineL[5],lineL[6]+"%"]))
			
			outline = "\t".join(outlineL)
			outL.append(outline)
	out = "\n".join(outL)
	print "buhege:" + str(NotHegeCount)
	
	return out



	


def main():
	iSNV_F = iSNVTableP + "iSNV_with_SNP.all.txt"
	iSNVinfoF = iSNVTableP + "iSNV_info.txt"
	SNPinfoF = iSNVTableP + "SNP_info.txt"	

	SAMPIDs = sample_ID(iSNV_F)
	print SAMPIDs


	for SAMPID in SAMPIDs:
		print SAMPID
		outvcfF = outvcfP + "/" + SAMPID + ".iSNV.vcf"
		Fexit = os.path.exists(outvcfF)
		if Fexit == True:
			os.remove(outvcfF)
		header = vcf_Header(SAMPID)
		outvcfFO = open(outvcfF,'a')
		outvcfFO.write(header + "\n")
		outvcfFO.close()
	
		

		sampPosL = iSNVsampPos(iSNV_F,SAMPID)		
		print "sampPosL" + "\t" + str(len(sampPosL))
			
		snvInfoD = samp_posi_info(iSNVinfoF,sampPosL,SAMPID)
		SNPInfoD = samp_posi_info(SNPinfoF,sampPosL,SAMPID)

	


		out = vcfFomat(sampPosL,snvInfoD,SNPInfoD)
		outvcfFO = open(outvcfF,'a')
		outvcfFO.write(out + "\n")
		outvcfFO.close()



if __name__ == "__main__":
	#sample = int(sys.argv[1])

	#iSNVTableP = "/home/liuhj/liuhj_2_44/test/table/"
	#outvcfF = "/home/liuhj/liuhj_2_44/test/" + str(sample) + ".outF.vcf"
	#SAMPID = str(sample) + "_2_26695"
	#RefChromID = "ref_26695"

	iSNVTableP = str(sys.argv[1])
	outvcfP = str(sys.argv[2])	
	RefChromID = str(sys.argv[3])
	#SAMPID = str(sys.argv[3])

	FREQ_threshold = 1
	baseCovThreshold = 10
	#baseFreqThreds = 0.05
	

	main()












