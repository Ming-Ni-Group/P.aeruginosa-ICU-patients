
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
			
			HeaderL.append(HeaderID.split(".sort")[0])
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
			Pos = int(lL[0])
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
					if Freq >= 0.05:
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
				if FREQ >= minSNPFreq:
					count += 1
		if count >= float(minSamppercent) * len(samples):
			personCommonPosi[posi] = ''
			
	return personCommonPosi



def readInfo(InfoF):
	cladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[1]
			clade = "-".join(InfoFl_Tags[2].split("_"))
			if clade  not in cladeSampDic:
				cladeSampDic[clade] = []
			cladeSampDic[clade].append(seqID)				
	return cladeSampDic




def HomoSites(HomoSitesF):
	HomoSitesDic = {}
	for HomoSitesFl in open(HomoSitesF).readlines():
		if HomoSitesFl != "\n":
			HomoSitesDic[HomoSitesFl.split("\n")[0].split("\t")[0]] = ''
	return HomoSitesDic


def RefRepeatMaskerRegions(RepeatMaskedGnm):
	for RepeatMaskedGnml in open(RepeatMaskedGnm).readlines():
		if ">" in RepeatMaskedGnml:
			gnm = ''
		else:
			gnm += RepeatMaskedGnml.split("\n")[0]
	repeatRegPosi={}
	for gnmPosiIndex in range(0,len(gnm)):
		if gnm[gnmPosiIndex] == "N":
			repeatRegPosi[gnmPosiIndex+1] = ''
	return repeatRegPosi


def FiltRefRepeatRegions(RepeatPosiF):
	repeatRegPosi = {}
	for RepeatPosiFl in open(RepeatPosiF).readlines():
		repeatRegPosi[int(RepeatPosiFl.split("\t")[0])] = ''
	return repeatRegPosi



def main():
	cladeSampDic = readInfo(InfoF)
	HomoSitesDic = HomoSites(HomoSitesF)
	repeatRegPosi = FiltRefRepeatRegions(RepeatF)
	print("clade samps :")
	print(cladeSampDic)
	print("Homo Sites num :" + str(len(HomoSitesDic.keys())))
	print("Repeat Sites num :" + str(len(repeatRegPosi.keys())))

	cladeLst = list(cladeSampDic.keys())
	cladeLst.sort()
	CladeCommSites_Dic = {}
	for Clade in cladeSampDic:
		iSNV_SNP_table = variants_P + "/" + Clade + "/all.iSNV_with_SNP.pyResults.txt"
		lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
		print("Clade: " + Clade)
		print(lieD)
		cladeSamps = cladeSampDic[Clade]
		print(cladeSamps)

		samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(cladeSamps,lieD,lsD)
		personNAPosi = Person_NA_Posi(cladeSamps,sampNA_PosiDic)
		personCommonPosi = Person_commonSNP_Posi(cladeSamps,samp_FreqDic,lsD)
		CladeCommSites_Dic[Clade] = personCommonPosi

	countDic = {}
	for Clade1 in cladeLst:
		for Clade2 in cladeLst:
			clade1SNPSites = CladeCommSites_Dic[Clade1]
			clade2SNPSites = CladeCommSites_Dic[Clade2]
			comm,diff = 0,0
			for Site1 in clade1SNPSites:
				if Site1 in clade2SNPSites:
					comm += 1
				else:
					diff += 1
			countDic[Clade1+"--"+Clade2] = {"comm":comm,"diff":diff}

	commOut,diffOut = [],[]
	header = ["clade"]
	for Clade1 in cladeLst:
		header.append(Clade1)
		commOutl,diffOutl = [Clade1],[Clade1]
		for Clade2 in cladeLst:
			commOutl.append(str(countDic[Clade1+"--"+Clade2]["comm"]))
			diffOutl.append(str(countDic[Clade1+"--"+Clade2]["diff"]))
		commOut.append("\t".join(commOutl))
		diffOut.append("\t".join(diffOutl))

	out_commF = out_P + "/SharedSNP_betweenDiffClades.txt"
	if ( os.path.exists(out_commF)):
		os.remove(out_commF)
	out_commF_O=open(out_commF,'a')
	out_commF_O.write("\t".join(header) + "\n")
	out_commF_O.write("\n".join(commOut))
	out_commF_O.close()	

	out_diffF = out_P + "/DiffSNP_betweenDiffClades.txt"
	if ( os.path.exists(out_diffF)):
		os.remove(out_diffF)
	out_commF_O=open(out_diffF,'a')
	out_commF_O.write("\t".join(header) + "\n")
	out_commF_O.write("\n".join(diffOut))
	out_commF_O.close()	


	repeatRegPosi = FiltRefRepeatRegions(RepeatF)
	RepeatMaskerSites = RefRepeatMaskerRegions(RepeatMaskerF)
	B,M = 0,0
	for RepeatMaskerSite in RepeatMaskerSites:
		if RepeatMaskerSite in repeatRegPosi:
			B += 1
		else:
			M += 1
	R = len(repeatRegPosi) - B

	print("\t".join(["all_R","all_M","single_R","single_M","Both"]))
	print("\t".join([str(len(repeatRegPosi)),str(len(RepeatMaskerSites)),str(R),str(M),str(B)]))


if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--variants_P",
					  dest = "variants_P",
					  default = "",
					  metavar = "path",
					  help = "iSNV table path.  [required]")
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

	parser.add_option("-R","--RepeatF",
					  dest = "RepeatF",
					  default = "",
					  metavar = "file",
					  help = "Repeat Sites file.  [required]")
	parser.add_option("-G","--RepeatMaskerF",
					  dest = "RepeatMaskerF",
					  default = "",
					  metavar = "file",
					  help = "RepeatMasker genome file.  [required]")


	parser.add_option("-M","--minSNPFreq",
					  dest = "minSNPFreq",
					  default = "0.9",
					  metavar = "float",
					  help = "min Freq of SNP (0-100%).  [required]")
	parser.add_option("-S","--minSamppercent",
					  dest = "minSamppercent",
					  default = "1",
					  metavar = "float",
					  help = "min percent of Samples of SNP (0-100%).  [required]")

	parser.add_option("-H","--HomoSitesF",
					  dest = "HomoSitesF",
					  default = "",
					  metavar = "file",
					  help = "Filt Ecoli homo Posi file.  [required]")


	(options,args) = parser.parse_args()
	variants_P   = os.path.abspath(options.variants_P)
	out_P	= os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	RepeatF = os.path.abspath(options.RepeatF)
	HomoSitesF   = os.path.abspath(options.HomoSitesF)
	RepeatMaskerF = os.path.abspath(options.RepeatMaskerF)

	minSNPFreq	 = float(options.minSNPFreq)
	minSamppercent	 = float(options.minSamppercent)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




