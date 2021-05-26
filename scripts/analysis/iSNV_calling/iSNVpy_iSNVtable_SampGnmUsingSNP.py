
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
	NADic = {}

	col = lieD[sample]
	
	for Posi in lsD:
		Freq = lsD[Posi][col]
		if Freq != "NA" :
			if Freq != "NO":				
				Freq = float(lsD[Posi][col])
				if Freq >= MinFreq:
					#print(Freq)
					samp_FreqDic[Posi] = Freq
		else:
			NADic[Posi]=''
			
	return samp_FreqDic,NADic





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





def reads_AltBaseDic(sample,samp_FreqDic,ntfreqlieD,ntfreqlsD):
	outLst = []
	sampSiteInfoD = {}
	SampCol = ntfreqlieD[sample]
	AlleDic = {}
	for Cite in samp_FreqDic:
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
		AlleDic[Cite] = ALT
	return AlleDic


	



def subSNPGnm(sample,samp_FreqDic,NADic,refID,RefGnm,AlleDic):
	RefChromID = refID.split(">")[1]

	Gnm = ''
	Ncount = 0
	for CiteIndex in range(0,len(RefGnm)):
		Cite = CiteIndex + 1
		if Cite in samp_FreqDic:
			Base = AlleDic[Cite]
		elif Cite in NADic:
			Base = "N"
			Ncount += 1
		else:
			Base = RefGnm[CiteIndex]

		Gnm += Base
	return Gnm,Ncount







def main():

	refID,RefGnm = readRefRefGnm(RefRefGnmF) 

	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	ntfreqlieD,ntfreqlsD = Read_ntfreqTable(ntfreqTable)
	

	sampLst = []
	for samp in lieD:
		if samp != "#Posi" and samp != "snv/SNP":
			sampLst.append(samp)

	person = iSNV_SNP_table.split("/")[-2]

	GnmCountF = out_P + "/" + person + ".subSNPGnm.NCount.fasta"
	if ( os.path.exists(GnmCountF)):
		os.remove(GnmCountF)


	for sample in sampLst:
		print("caculating for sample : " + sample)
		samp_FreqDic,NADic = samp_Freq_Dic(sample,lieD,lsD,MinFreq)
		#print(samp_FreqDic)
		AlleDic = reads_AltBaseDic(sample,samp_FreqDic,ntfreqlieD,ntfreqlsD)
		#print(AlleDic)
		Gnm,Ncount = subSNPGnm(sample,samp_FreqDic,NADic,refID,RefGnm,AlleDic)


		GnmF = out_P + "/" + sample.split("_BDM")[0].split("BJ13-")[1] + ".minFreq" + str(MinFreq) + ".subSNPGnm.fasta"
		if ( os.path.exists(GnmF)):
			os.remove(GnmF)
		out_snv_F_O=open(GnmF,'a')
		out_snv_F_O.write(">" + sample.split("_BDM")[0].split("BJ13-")[1] + "\n")
		out_snv_F_O.write(Gnm + "\n")
		#out_snv_F_O.write(refID + "\n")
		#out_snv_F_O.write(RefGnm + "\n")
		out_snv_F_O.close()	

		print(sample + "       ----  Done. ")



		out_snv_F_O=open(GnmCountF,'a')
		out_snv_F_O.write(sample.split("_BDM")[0].split("BJ13-")[1] +"\t" + str(Ncount) + "\n")
		out_snv_F_O.close()	





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
					  default = "0.5",
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




