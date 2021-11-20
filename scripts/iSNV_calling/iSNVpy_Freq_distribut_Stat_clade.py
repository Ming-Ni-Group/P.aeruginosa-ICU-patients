
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def iSNV_SNP_tableRead(iSNV_SNP_table):
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
	HeaderLs = iSNVSNPtablels[0].split("\n")[0].split("\t")
	HeaderL = []
	for HeaderID in HeaderLs:
		if HeaderID != "#Posi" and HeaderID !=  "snv/SNP":
			
			HeaderL.append(HeaderID.split("_BDM")[0].split("BJ13-")[1])
		else:
			HeaderL.append(HeaderID)
	lie = 0
	lieD = {}
	for Header in HeaderL:
		lieD[Header] = lie
		lie += 1
	lsD = {}
	for iSNVSNPtablel in iSNVSNPtablels:
		if "#" not in iSNVSNPtablel and iSNVSNPtablel != "\n" :
			lL = iSNVSNPtablel.split("\n")[0].split("\t")
			Pos = lL[0]
			lsD[Pos] = lL
	return lieD,lsD



def samp_Freq_Dic(samples,lieD,lsD):
	samp_FreqDic = {}
	sampNA_PosiDic = {}
	print(lieD)
	for sample in samples:
		Lie = lieD[sample]
		samp_FreqDic[sample] = {}
		sampNA_PosiDic[sample] = []
		for Posi in lsD:
			Freq = lsD[Posi][Lie]
			if Freq != "NA" :
				if Freq != "NO":				
					Freq = float(lsD[Posi][Lie])
					if Freq >= 0.02:
						samp_FreqDic[sample][Posi] = Freq
				else:
					Freq = 0
			else:
				sampNA_PosiDic[sample].append(Posi)
	#print(samp_FreqDic.keys())
	return samp_FreqDic,sampNA_PosiDic



def Person_NA_Posi(samples,sampNA_PosiDic):
	personNAPosi = {}
	for sample in samples:
		sampleNAPosi = sampNA_PosiDic[sample]
		for posi in sampleNAPosi:
			if posi not in personNAPosi:
				personNAPosi[posi] = ''
	return personNAPosi


def Person_commonSNP_Posi(samples,samp_FreqDic,lsD):
	personCommonPosi = {}
	
	for posi in lsD.keys():
		count = 0
		for sample in samples:
			sampleFreqDic = samp_FreqDic[sample]
			if posi in sampleFreqDic:
				FREQ = sampleFreqDic[posi]
				if FREQ >= snv_MaxFreq:
					count += 1
		if count == len(samples):
			personCommonPosi[posi] = ''
			
	return personCommonPosi




def readInfo(InfoF,person):
	print(person)
	cladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[0]
			clade = InfoFl_Tags[20]

			if person in clade:
				if clade  not in cladeSampDic:
					cladeSampDic[clade] = []
				cladeSampDic[clade].append(seqID)
	return cladeSampDic





def tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,snvSNPcount_Dic,personSNP_Lst,personsnv_Lst):

	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}
	#print(len(personCommonPosi))
	#print("len common")
	#for XX in range(0,101,step):
		#tongji_Dic[XX] = 0

	snvCount = 0
	SNPCount = 0

	
	for Posi in FreqDic:
		if Posi not in personNAPosi and Posi not in personCommonPosi:
			Freq = FreqDic[Posi]

			if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
				snvCount += 1
				if Posi not in personsnv_Lst:
					personsnv_Lst.append(Posi)

			elif Freq >= snv_MaxFreq :
				SNPCount += 1
				if Posi not in personSNP_Lst:
					personSNP_Lst.append(Posi)


	print("\t".join([sample,str(snvCount),str(SNPCount)]))
	snvSNPcount_Dic[sample] = "\t".join([sample,str(snvCount),str(SNPCount),str(len(personCommonPosi))])

	return snvSNPcount_Dic,personSNP_Lst,personsnv_Lst






'''
	freq_flagLst=[]
	snvCount_list=[]
	outlst = []
	for XX in range(0,101,step):
		outlst.append(sample +"\t" + str(XX) + "\t" + str(tongji_Dic[XX]))
		freq_flagLst.append(XX)
		snvCount_list.append(tongji_Dic[XX])

	outlines = "\n".join(outlst) 
	out_F_O=open(out_F,'a')
	out_F_O.write(outlines)
	out_F_O.close()	
	
	out_plotF = out_F.split(".txt")[0] + ".pdf"
	plt.figure(dpi=300,figsize=(8,3))
	plt.bar(range(len(snvCount_list)), snvCount_list,color='red',tick_label=freq_flagLst)  


	plt.title(sample.split(".sort")[0])
	plt.xlabel('Freq', fontsize = 13)
	plt.ylabel('count', fontsize = 13)
	plt.xticks(rotation=90)
	plt.savefig(out_plotF)
	plt.savefig(out_plotF.split(".pdf")[0] + ".png")
'''

## snv and SNP count 
	#out_snvSNPCountF = out_F.split(".txt")[0] + ".snvSNPcount.txt"
	#out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	#out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount"]) + "\n")
	#out_snvSNPCountF_O.write("\t".join([sample,str(snvCount),str(SNPCount)]))
	#out_snvSNPCountF_O.close()	

	

def Posi_lines(cladeSamps,lieD,lsD,Cite_Lst):
	SampLst = []
	for cladeSamp in cladeSamps:
		SampLst.append("S" + cladeSamp)
	outlines = ["\t".join(["#posi","no.SNP"] + SampLst)]
	for POSI in Cite_Lst:
		lineLst = [POSI,'']
		SNPno = 0
		for cladeSamp in cladeSamps:
			Samp_col = lieD[cladeSamp]
			Freq = lsD[POSI][Samp_col]
			lineLst.append(Freq)
			if Freq != "NO":
				SNPno += 1
		lineLst[1] = str(SNPno)
		line = "\t".join(lineLst)
		outlines.append(line)

	lines = "\n".join(outlines)
	return lines









def main():
	person = iSNV_SNP_table.split("/")[-2]
	if person == "P5_withP2":
		person = "P5"
	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)


	#iSNVtable_lst = list(lieD.keys())
	cladeSampDic = readInfo(InfoF,person)
	print(person)
	print(cladeSampDic)
	for clade in cladeSampDic:
		cladeSamps = cladeSampDic[clade]
		print(cladeSamps)


		samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(cladeSamps,lieD,lsD)
		personNAPosi = Person_NA_Posi(cladeSamps,sampNA_PosiDic)
		personCommonPosi = Person_commonSNP_Posi(cladeSamps,samp_FreqDic,lsD)
		snvSNPcount_Dic = {}
		personSNP_Lst = []	
		personsnv_Lst = []	
		for sample in cladeSamps:
			#print(sample)

			snvSNPcount_Dic,personSNP_Lst,personsnv_Lst = tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,snvSNPcount_Dic,personSNP_Lst,personsnv_Lst)		
		#person = iSNV_SNP_table.split("/")[-2]
		out_snvSNPCountF = out_P + "/" + clade + ".snvSNPcount.txt"
		out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
		out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount","personCommonCitesNum"]) + "\n")
		for sampleID in snvSNPcount_Dic:	
			out_snvSNPCountF_O.write(snvSNPcount_Dic[sampleID] + "\n")
		out_snvSNPCountF_O.close()	

		personSNP_Lst.sort()
		SNPCiteslines = Posi_lines(cladeSamps,lieD,lsD,personSNP_Lst)
		out_SNP_F = out_P + "/" + clade + ".SNP_Posi_list.txt"
		if ( os.path.exists(out_SNP_F)):
			os.remove(out_SNP_F)
		out_SNP_F_O=open(out_SNP_F,'a')
		out_SNP_F_O.write(SNPCiteslines + "\n")
		out_SNP_F_O.close()	

		print(len(personSNP_Lst))
		print(len(personsnv_Lst))
		personsnv_Lst.sort()
		snvCiteslines = Posi_lines(cladeSamps,lieD,lsD,personsnv_Lst)
		out_snv_F = out_P + "/" + clade + ".snv_Posi_list.txt"
		if ( os.path.exists(out_snv_F)):
			os.remove(out_snv_F)
		out_snv_F_O=open(out_snv_F,'a')
		out_snv_F_O.write(snvCiteslines + "\n")
		out_snv_F_O.close()	








if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table",
					  dest = "iSNV_SNP_table",
					  default = "",
					  metavar = "file",
					  help = "iSNV table.  [required]")
	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-f","--InfoF",
					  dest = "InfoF",
					  default = "",
					  metavar = "file",
					  help = "Info File.  [required]")


	parser.add_option("-m","--snv_MinFreq",
					  dest = "snv_MinFreq",
					  default = "0.05",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--snv_MaxFreq",
					  dest = "snv_MaxFreq",
					  default = "0.95",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")


	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	snv_MinFreq		= float(options.snv_MinFreq)
	snv_MaxFreq		= float(options.snv_MaxFreq)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




