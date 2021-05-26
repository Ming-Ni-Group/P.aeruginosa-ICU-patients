
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


def CovtableRead(Covtable):
	Covtablels = open(Covtable,'r').readlines()
	HeaderLs = Covtablels[0].split("\n")[0].split("\t")
	HeaderL = []
	for HeaderID in HeaderLs:
		HeaderL.append(HeaderID)
	lie = 0
	CvlieD = {}
	for Header in HeaderL:
		CvlieD[Header] = lie
		lie += 1
	CvlsD = {}
	for l in Covtablels:
		if "#" not in l and l != "\n" :
			lL = l.split("\n")[0].split("\t")
			Pos = int(lL[0])
			CvlsD[Pos] = lL
	return CvlieD,CvlsD



def FiltPosi_HomowithOtherBacteria(FiltOtherBacteriaHomoPosiF):
	EcoliFiltedPosiDic = {}
	for FiltOtherBacteriaHomoPosiFl in open(FiltOtherBacteriaHomoPosiF).readlines():
		if FiltOtherBacteriaHomoPosiFl != "\n":
			EcoliFiltedPosiDic[int(FiltOtherBacteriaHomoPosiFl.split("\n")[0].split("\t")[0])] = ''
	return EcoliFiltedPosiDic


def FiltRefRepeatRegions(RepeatPosiF):
	repeatRegPosi = {}
	for RepeatPosiFl in open(RepeatPosiF).readlines():
		repeatRegPosi[int(RepeatPosiFl.split("\t")[0])] = ''
	return repeatRegPosi



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



def readInfo(InfoF,person):
	print(person)
	cladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[18]
			clade = InfoFl_Tags[20]
			#if seqID not in ["P1-S1","P1-S2","P1-S3","P2-S1","P5-S17","P7-S3"]:	
			if person in clade:
				if clade  not in cladeSampDic:
					cladeSampDic[clade] = []
				cladeSampDic[clade].append(seqID)			
		
	return cladeSampDic



def readSnpEffF(sample,sampleSnpEffP,RefDic,AltDic):
	AltDic[sample] = {}	
	sampleSnpEffF = sampleSnpEffP + "/" + sample + ".sort.rmdup.mpileup.minFreq0.05.snpeff.vcf"
	for sampleSnpEffFl in open(sampleSnpEffF).readlines():
		if "#" not in sampleSnpEffFl and sampleSnpEffFl != "\n":
			sampleSnpEffFl_Tags = sampleSnpEffFl.split("\n")[0].split("\t")
			Posi = int(sampleSnpEffFl_Tags[1])			
			ref = sampleSnpEffFl_Tags[3]
			alt = sampleSnpEffFl_Tags[4]
			if Posi not in RefDic:
				RefDic[Posi] = ref
			AltDic[sample][Posi] = alt

	return RefDic,AltDic




def main():
	cladesLst = ["P1","P3-1","P3-2","P3-3","P4","P5-1","P5-1","P5-1","P6","P7"]
	repeatRegPosi = FiltRefRepeatRegions(RepeatPosiF)
	HomoPosiOtherDic = FiltPosi_HomowithOtherBacteria(FiltOtherBacteriaHomoPosiF)

	AllSamps = []
	AllPosi,allNAPosi = {},{}
	RefDic,AltDic = {},{}

	for clade in cladesLst:
		iSNVtable = iSNV_SNP_table_P + "/" + clade + "/" + "all.iSNV_with_SNP.pyResults.txt"
		Covtable = iSNV_SNP_table_P + "/" + clade + "/" + "out.samps-Coverage.bigtable.txt"
		lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
		CvlieD,CvlsD = CovtableRead(Covtable)
		for posi in lsD:
			if posi not in AllPosi:
				AllPosi[posi] = ''
		for CvPos in CvlsD:
			for CvSamp in CvlieD:
				CoSampLie = CvlieD[CvSamp]
				CovSum = int(CvlsD[CvPos][CoSampLie])
				if CovSum < 50 and CvPos not in allNAPosi:
					allNAPosi[CvPos] = ''

	AllDepreliPos = {}
	for Posi in AllPosi:
		if Posi not in allNAPosi:
			AllDepreliPos[Posi] = ''
	Depreli_OtherHom,Depreli_repeat = {},{}
	finalPosi = {}
	for DepreliPos in AllDepreliPos:
		if DepreliPos in HomoPosiOtherDic:
			Depreli_OtherHom[DepreliPos] = ''
		elif DepreliPos in repeatRegPosi:
			Depreli_repeat[DepreliPos] = ''
		if posi not in allNAPosi and posi not in Depreli_OtherHom and posi not in Depreli_repeat:
			finalPosi[posi] = ''
	
	RefDic,AltDic = {},{}
	for clade in cladesLst:
		iSNVtable = iSNV_SNP_table_P + "/" + clade + "/" + "all.iSNV_with_SNP.pyResults.txt"
		lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
		for posi in finalPosi:
			for samp in lieD:
				readSnpEffF(sample,sampleSnpEffP,RefDic,AltDic)	
				if samp not in AllSamps:
					AllSamps.append(samp)


	SampSNPPosi = {}
	for clade in cladesLst:
		iSNVtable = iSNV_SNP_table_P + "/" + clade + "/" + "all.iSNV_with_SNP.pyResults.txt"
		lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)
		for posi in finalPosi:
			for samp in lieD:
				SampSNPPosi[samp] = {}
				Freq = float(lsD[posi][lieD[samp]])
				if Freq >= SNP_minFreq:
					SampSNPPosi[samp][posi] = AltDic[posi]

	print("allSNP")
	outs = []
	finalPosi.sort()
	for samp in AllSamps:
		sampSNPSeq = ''
		for posi in finalPosi:
			if posi in SampSNPPosi[samp]:
				base = SampSNPPosi[samp][posi] 
			else:
				base = RefDic[posi]
			sampSNPSeq += base 
		outs.append(sampSNPSeq)
	out = "\n".join(outs)
				
	out_SNP_F = out_P + "/samps.SNPseq.-OtherBac-repeat.fasta"
	if ( os.path.exists(out_SNP_F)):
		os.remove(out_SNP_F)
	out_SNP_F_O=open(out_SNP_F,'a')
	out_SNP_F_O.write(out + "\n")
	out_SNP_F_O.close()	




	header = "\t".join(["samp"] + AllSamps)
	ls = []
	for samp1 in AllSamps:
		l = [samp1]
		for samp2 in AllSamps:
			diffC,commC=0,0
			for posi in finalPosi:
				if posi in SampSNPPosi[samp1]:
					base1 = SampSNPPosi[samp1][posi] 
				else:
					base1 = ''
				
				if posi in SampSNPPosi[samp2]:
					base2 = SampSNPPosi[samp2][posi] 
				else:
					base2 = ''

				if base1 != '' or base1 != '':
					if (base1 != base2) :
						diffC += 1
					if (base1 != base2) :
						commC += 1
			l.append(str(diffC) + "/" + str(commC))
		ls.append("\t".join(l))
	out = "\n".join(ls)

	out_F = out_P + "/samps.diff-common.Count.txt"
	if ( os.path.exists(out_F)):
		os.remove(out_F)
	out_SNP_F_O=open(out_F,'a')
	out_SNP_F_O.write(out + "\n")
	out_SNP_F_O.close()	
			







if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--iSNV_SNP_table_P",
					  dest = "iSNV_SNP_table_P",
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
	parser.add_option("-s","--sampleSnpEffP",
					  dest = "sampleSnpEffP",
					  default = "",
					  metavar = "path",
					  help = "sample SnpEff Path.  [required]")

	parser.add_option("-R","--RepeatRegion",
					  dest = "RepeatRegion",
					  default = "",
					  metavar = "file",
					  help = "Repeat Regions .  [required]")
	parser.add_option("-E","--FiltOtherBacteriaHomoPosiF",
					  dest = "FiltOtherBacteriaHomoPosiF",
					  default = "",
					  metavar = "file",
					  help = "Filt Ecoli homo Posi file.  [required]")


	parser.add_option("-M","--SNP_minFreq",
					  dest = "SNP_minFreq",
					  default = "0.9",
					  metavar = "float",
					  help = "min Freq of SNP to generate the sequences to build the tree (0-100%).  [required]")



	(options,args) = parser.parse_args()
	iSNV_SNP_table_P   = os.path.abspath(options.iSNV_SNP_table_P)
	out_P		  = os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	sampleSnpEffP = os.path.abspath(options.sampleSnpEffP)
	RepeatMaskedGnm = os.path.abspath(options.RepeatMaskedGnm)
	FiltOtherBacteriaHomoPosiF   = os.path.abspath(options.FiltOtherBacteriaHomoPosiF)

	SNP_minFreq		= float(options.SNP_minFreq)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




