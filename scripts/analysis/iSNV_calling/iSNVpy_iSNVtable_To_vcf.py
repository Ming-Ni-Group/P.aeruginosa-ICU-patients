
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def readRefRefGnm(RefRefGnmF):
	print("Loading Ref genome file :" + RefRefGnmF)
	RefGnm = ''
	for RefRefGnmFl in open(RefRefGnmF).readlines():
		if ">" in RefRefGnmFl:
			refID = RefRefGnmFl.split("\n")[0]
		else:
			RefGnm += RefRefGnmFl.split("\n")[0]
	print("ref genome loaded!" )
	print("ref genome base num: " + str(len(RefGnm)))
	print()
	return refID,RefGnm



def iSNV_SNP_tableRead(iSNV_SNP_table):
	print("Loading iSNV table file :" + iSNV_SNP_table)
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
	HeaderLs = iSNVSNPtablels[0].split("\n")[0].split("\t")
	lie = 0
	lieD = {}
	for Header in HeaderLs:
		lieD[Header] = lie
		lie += 1
	lsD = {}
	for iSNVSNPtablel in iSNVSNPtablels:
		if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
			lL = iSNVSNPtablel.split("\n")[0].split("\t")
			Pos = int(lL[0])
			lsD[Pos] = lL
	print("iSNV table loaded!" )
	print("iSNV table samples num: " + str(len(lieD.keys()) - 2 ))
	print()
	return lieD,lsD


def samp_Freq_Dic(sample,lieD,lsD,MinFreq):
	samp_FreqDic = {}
	col = lieD[sample]
	for Posi in lsD:
		Freq = lsD[Posi][col]
		if Freq != "NA" :
			if Freq != "NO":				
				Freq = float(lsD[Posi][col])
				if Freq >= MinFreq:
					#print(Freq)
					samp_FreqDic[Posi] = Freq
			
	return samp_FreqDic






def vcf_Header(sample):
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
	#HeaderL.append('##FORMAT=<ID=MDP,Number=1,Type=Integer,Description="Raw Read Depth from mpileup">')
	HeaderL.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= MapQThreds">')
	HeaderL.append('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">')
	HeaderL.append('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">')
	HeaderL.append('##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">')
	HeaderL.append('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	' + sample)

	HEADER = "\n".join(HeaderL)
	return HEADER




def Read_ntfreqTable(ntfreqTable):
	print("Loading ntfreq big Table file :" + ntfreqTable)

	ntfreqTablels = open(ntfreqTable,'r').readlines()
	HeaderLs = ntfreqTablels[0].split("\n")[0].split("\t")
	lie = 0
	ntfreqlieD = {}
	for Header in HeaderLs:
		ntfreqlieD[Header] = lie
		lie += 1
	ntfreqlsD = {}
	for ntfreqTablel in ntfreqTablels:
		if "#" not in ntfreqTablel and ntfreqTablel != "\n" :
			lL = ntfreqTablel.split("\n")[0].split("\t")
			Pos = int(lL[0])
			ntfreqlsD[Pos] = lL
	print("ntfreq big Table loaded!" )
	print()
	return ntfreqlieD,ntfreqlsD









def convertToVcf(sample,samp_FreqDic,ntfreqlieD,ntfreqlsD,refID,RefGnm):
	outLst = []
	sampSiteInfoD = {}
	SampCol = ntfreqlieD[sample]
	RefChromID = refID.split(">")[1]

	for Cite in samp_FreqDic:
		REF = RefGnm[Cite-1]
		ALT = ''
		Samp_ntfreq = ntfreqlsD[Cite][SampCol]
		nt_pattern = Samp_ntfreq.split("_")
		Anum=int(nt_pattern[0])
		Gnum=int(nt_pattern[1])
		Cnum=int(nt_pattern[2])
		Tnum=int(nt_pattern[3])
		totnum=int(nt_pattern[4])

		charsCount = 0
		NumD = {"A":Anum,"G":Gnum,"C":Cnum,"T":Tnum}
		
		FreqD = {}

		Afreq = round(float(Anum)/totnum,4)
		FreqD["A"] = Afreq
		Gfreq = round(float(Gnum)/totnum,4)
		FreqD["G"] = Gfreq
		Cfreq = round(float(Cnum)/totnum,4)
		FreqD["C"] = Cfreq
		Tfreq = round(float(Tnum)/totnum,4)
		FreqD["T"] = Tfreq

		for Base in FreqD:
			if FreqD[Base] == samp_FreqDic[Cite] :
				ALT = Base

	
		RD = str(NumD[REF])
		AD = str(NumD[ALT])
		FREQ = samp_FreqDic[Cite]
		DP = totnum
		if len(str(FREQ)) > 5:
			print(FREQ)

		outlineL = []
		outlineL.append(RefChromID)
		outlineL.append(str(Cite))
		outlineL.append(".")
		outlineL.append(REF)  ##REF
		outlineL.append(ALT)  ##ALT
		outlineL.append(".")
		outlineL.append("PASS")
		outlineL.append(";".join(["ADP="+ str(DP),"WT=0","HET=0","HOM=1","NC=0"]))
		outlineL.append("GT:GQ:DP:RD:AD:FREQ")
		outlineL.append(":".join(["1/1","255",str(DP),str(RD),str(AD),str(FREQ*100)+"%"]))			
		outline = "\t".join(outlineL)

		outLst.append(outline)
	out = "\n".join(outLst)


	return out				




def main():

	refID,RefGnm = readRefRefGnm(RefRefGnmF) 

	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	ntfreqlieD,ntfreqlsD = Read_ntfreqTable(ntfreqTable)
	

	sampLst = []
	for samp in lieD:
		if samp != "#Posi" and samp != "snv/SNP":
			sampLst.append(samp)
	for sample in sampLst:
		
		print("caculating for sample : " + sample)
		samp_FreqDic = samp_Freq_Dic(sample,lieD,lsD,MinFreq)
		vcflines = convertToVcf(sample,samp_FreqDic,ntfreqlieD,ntfreqlsD,refID,RefGnm)

		vcfHead = vcf_Header(sample)
		out_vcfF = out_P + "/" + sample + ".minFreq" + str(MinFreq) + ".iSNVpy.vcf"
		if ( os.path.exists(out_vcfF)):
			os.remove(out_vcfF)
		out_snv_F_O=open(out_vcfF,'a')
		out_snv_F_O.write(vcfHead + "\n")
		out_snv_F_O.write(vcflines + "\n")
		out_snv_F_O.close()	
		print(sample + "       ----  Done. ")








if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table",
					  dest = "iSNV_SNP_table",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-n","--ntfreqTable",
					  dest = "ntfreqTable",
					  default = "",
					  metavar = "file",
					  help = "ntfreq big table  [required]")
	parser.add_option("-r","--RefRefGnmF",
					  dest = "RefRefGnmF",
					  default = "",
					  metavar = "file",
					  help = "ref genome file .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")




	parser.add_option("-m","--MinFreq",
					  dest = "MinFreq",
					  default = "0.05",
					  metavar = "float",
					  help = "min Freq of cite to consit in Vcf (0-100%).  [optional]")





	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	ntfreqTable   = os.path.abspath(options.ntfreqTable)
	RefRefGnmF   = os.path.abspath(options.RefRefGnmF)
	out_P		  = os.path.abspath(options.out_P)
	MinFreq		= float(options.MinFreq)


	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




