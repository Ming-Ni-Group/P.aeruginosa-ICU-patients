
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 


def iSNV_SNP_tableRead(iSNV_SNP_table):
	iSNVSNPtablels = open(iSNV_SNP_table,'r').readlines()
	HeaderL = iSNVSNPtablels[0].split("\n")[0].split("\t")
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
	#print(lieD)
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
				if FREQ >= 0.95:
					count += 1
		if count == len(samples):
			personCommonPosi[posi] = ''
			
	return personCommonPosi




def tongjiFREQ(sample,samp_FreqDic,out_F,personNAPosi,personCommonPosi,snvSNPcount_Dic):

	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}
	#print(len(personCommonPosi))
	#print("len common")
	for XX in range(0,101,step):
		tongji_Dic[XX] = 0

	snvCount = 0
	SNPCount = 0

	for Posi in FreqDic:
		if Posi not in personNAPosi and Posi not in personCommonPosi:
			Freq = FreqDic[Posi]

			for XX in range(0,101,step):
				if Freq * 100 >= XX and Freq*100 < XX+step:
					tongji_Dic[XX]  += 1

			if Freq >= minFreq and Freq < 1-minFreq:
				snvCount += 1
			elif Freq >= 1-minFreq :
				SNPCount += 1

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


## snv and SNP count 
	#out_snvSNPCountF = out_F.split(".txt")[0] + ".snvSNPcount.txt"
	#out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	#out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount"]) + "\n")
	#out_snvSNPCountF_O.write("\t".join([sample,str(snvCount),str(SNPCount)]))
	#out_snvSNPCountF_O.close()	
	
	print("\t".join([sample,str(snvCount),str(SNPCount)]))
	snvSNPcount_Dic[sample] = "\t".join([sample,str(snvCount),str(SNPCount),str(len(personCommonPosi))])

	return snvSNPcount_Dic






def main():

	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
	#print(lieD)
	iSNVtable_lst = list(lieD.keys())
	sampIDs = []
	for sampID in iSNVtable_lst:
		if "#Posi" not in sampID and "snv/SNP" not in sampID:
			sampIDs.append(sampID)


	samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(sampIDs,lieD,lsD)
	personNAPosi = Person_NA_Posi(sampIDs,sampNA_PosiDic)
	personCommonPosi = Person_commonSNP_Posi(sampIDs,samp_FreqDic,lsD)
	snvSNPcount_Dic = {}

	for sample in sampIDs:
		#print(sample)
		out_F = out_P + "/" + sample + "iSNVpy.Freq_distribu.txt"
		snvSNPcount_Dic = tongjiFREQ(sample,samp_FreqDic,out_F,personNAPosi,personCommonPosi,snvSNPcount_Dic)

	person = iSNV_SNP_table.split("/")[-2]
	out_snvSNPCountF = out_P + "/" + person + ".snvSNPcount.txt"
	out_snvSNPCountF_O=open(out_snvSNPCountF,'a')
	out_snvSNPCountF_O.write("\t".join(["sample","snvCount","SNPCount","personCommonCitesNum"]) + "\n")
	for sampleID in snvSNPcount_Dic:	
		out_snvSNPCountF_O.write(snvSNPcount_Dic[sampleID] + "\n")
	out_snvSNPCountF_O.close()	





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

	parser.add_option("-m","--minFreq",
					  dest = "minFreq",
					  default = "0",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--maxFreq",
					  dest = "maxFreq",
					  default = "100",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-s","--step",
					  dest = "step",
					  default = "5",
					  metavar = "int",
					  help = "frame step to calculate.  [required]")


	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	minFreq		= float(options.minFreq)
	maxFreq		= float(options.maxFreq)
	step		   = int(options.step)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




