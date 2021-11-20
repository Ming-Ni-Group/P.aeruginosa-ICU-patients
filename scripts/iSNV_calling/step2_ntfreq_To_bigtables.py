
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def readntfreqFList_configFile(ntfreqFList_ntfreqFList_configF):
	ntfreqFList = []
	for configline in open(ntfreqFList_ntfreqFList_configF).readlines():
		if configline !="\n":
			ntfreqF = configline.split("\n")[0]
			ntfreqFList.append(ntfreqF)
	return ntfreqFList



def refgenome_Base(RefFile):
	print("loading reference genome: " + RefFile)
	refCiteBaseDic = {}
	seq = ''
	for RefFile in open(RefFile).readlines():
		if ">" in RefFile:
			refID = RefFile.split("\n")[0].split(">")[-1]
		if ">" not in RefFile:
			seq += RefFile.split("\n")[0]

	for Cite in range(0,len(seq)):
		posiIndex = Cite + 1
		refCiteBaseDic[posiIndex] = seq[Cite]
	print("ref genome loaded.")
	print("ref genome base Num : " + str(len(refCiteBaseDic)))
	print()
	return refID,refCiteBaseDic



def refgenome_Base11(RefFile):  ##may not be used
	print("loading reference genome: " + RefFile)
	refCiteBaseDic = {}
	refSeqDic = {}
	for RefFile in open(RefFile).readlines():
		if ">" in RefFile:
			refID = RefFile.split("\n")[0].split(">")[-1]
			refSeqDic[refID] = ''
		else:
			refSeqDic[refID] += RefFile.split("\n")[0]

	#print(refSeqDic)
	for refid in refSeqDic.keys():
		refCiteBaseDic[refid] = {}	
		for Cite in range(0,len(refSeqDic[refid])):
			posiIndex = Cite + 1
			refCiteBaseDic[refid][posiIndex] = refSeqDic[refid][Cite]
	print("ref genome loaded.")
	return refCiteBaseDic



def All_ntfreq_file(refID,ntfreqFLst,DEP_THRES):
	All_ntpattDic = {}
	All_StrandbiaDic = {}
	All_CovDic = {}
	All_ntfreqCitesLst = []

	file_num = 0
	Stat_Dic = {"IsSampleValid":{},"meanSiteDep":{},"validPosNo":{}}

	for ntfreqF in ntfreqFLst:
		file_num += 1
		print("loading ntfreq file " +  str(file_num) + "  :  " + ntfreqF)
		SampID = ntfreqF.split("/")[-1].split(".ntfreq")[0]
		All_ntpattDic[SampID] = {}
		All_StrandbiaDic[SampID] = {}
		All_CovDic[SampID] = {} 

		Valid_CiteNum = 0
		CovLst = []
		with open(ntfreqF, 'r') as file:
			for ntfreqFl in file:
				#print(ntfreqFl,)
		#for ntfreqFl in open(ntfreqF, 'r').readlines():				
			#AimTags = {}
				if "#" not in ntfreqFl:
					Tags = ntfreqFl.split("\t")
					refChrm = Tags[0]
					if refChrm != refID:
						print("ref genome ID not equel the refID in ntfreq file")
						break					
					Posi = int(Tags[1])
					#print(Posi)
					#if Posi not in All_ntfreqCitesLst:			
						#All_ntfreqCitesLst.append(Posi)
					#refbase = Tags[2]
					#if refbase != refCiteBaseDic[Posi]:
						#print("hhhhhhhhhhh")
					#ntfreqPatt = Tags[5]
					QCTotNum = int(Tags[10])
					#Anum = Tags[6]
					#Anum = Tags[7]
					#Anum = Tags[8]
					#Anum = Tags[9]
	
					#A_5_starnd = Tags[19]
					#G_5_starnd = Tags[20]
					#C_5_starnd = Tags[21]
					#T_5_starnd = Tags[22]
					
					CovLst.append(QCTotNum)
					All_CovDic[SampID][Posi] = QCTotNum
					if QCTotNum >= DEP_THRES:
						Valid_CiteNum += 1
					All_ntpattDic[SampID][Posi] = Tags[6]	+ "_" + Tags[7] + "_" + Tags[8] + "_" + Tags[9] + "_" + Tags[10]			
					All_StrandbiaDic[SampID][Posi] = Tags[19]+ "_" + Tags[20]+ "_" + Tags[21]+ "_" + Tags[22]
							

		if Valid_CiteNum >= VALID_SIZE_THRES:
			IsSampleValidFlag = 'yes'
		else:
			IsSampleValidFlag = 'no'

		Stat_Dic["IsSampleValid"][SampID] = IsSampleValidFlag
		Stat_Dic["meanSiteDep"][SampID] = np.mean(CovLst)
		Stat_Dic["validPosNo"][SampID] = Valid_CiteNum
	return All_ntpattDic,Stat_Dic,All_StrandbiaDic,All_CovDic,All_ntfreqCitesLst






def find_alle(refCiteBaseDic,All_ntpattDic,All_StrandbiaDic,DEP_THRES,MIN_DEP_THRES,STRANDED_RATIO_THRES,FREQ_THRES):
	AllAlle_Dic = {}
	AllAlleFreq_Dic = {}
	Strand_bia_index = {"A":0,"G":1,"C":2,"T":3}

	for Sample in All_ntpattDic.keys():
		AllAlle_Dic[Sample] = {}
		AllAlleFreq_Dic[Sample] = {}
		for Cite in All_ntpattDic[Sample]:
			#print(Cite)
			ntTags = All_ntpattDic[Sample][Cite].split("_")
			yesBaseLst = []
			yesBaseNumLst = []
			BaseLst = ["A","G","C","T"]
			BaseNumLst = [int(ntTags[0]),int(ntTags[1]),int(ntTags[2]),int(ntTags[3])]
			TotNum = int(ntTags[-1])

			countIndex = 0
			#print(BaseNumLst)
			for Base in BaseLst:		
				BaseNum = BaseNumLst[countIndex]
				#print(BaseNum)
				if BaseNum >= MIN_DEP_THRES:
					yesBaseLst.append(Base)
					yesBaseNumLst.append(BaseNum)
				
				countIndex += 1
			#print(Cite)
			#print(yesBaseNumLst)

			Flag = ''
			if yesBaseNumLst != []:				
				maxNum = max(yesBaseNumLst)
				maxIndex = yesBaseNumLst.index(max(yesBaseNumLst))  # 求最大值的索引
				maxBase = yesBaseLst[maxIndex]	
				
				Strand_bia = All_StrandbiaDic[Sample][Cite].split("_")[Strand_bia_index[maxBase]]
				Strand_bia_Flag = Strand_bia
				#Alle = "kkk"
				#AlleFreq = "kkk"
	
				if maxBase != refCiteBaseDic[Cite]:
					if TotNum >= DEP_THRES and maxNum >= MIN_DEP_THRES:
						if ( Strand_bia != "na" and (float(Strand_bia) >= STRANDED_RATIO_THRES ) and (float(Strand_bia) <= 1-STRANDED_RATIO_THRES)  ):
							Alle = maxBase
							AlleFreq = round(float(maxNum)/TotNum,4)
							if AlleFreq < float(FREQ_THRES):
								AlleFreq = "NO"
								Alle = "NO"
						else:
							Alle = "NO"
							AlleFreq = "NO"
							Flag = "max  depth 合格  链偏性不合格"
						#	print("ref:" + refCiteBaseDic[Cite])
						#	print(yesBaseLst)
						#	print(yesBaseNumLst)
						#	print(maxBase)
						#	print(maxNum)
						#	print(Alle)
						#	print(AlleFreq)
						#	print(Strand_bia)
						#	print()
	
					else:
						Alle = "NA"
						AlleFreq = "NA"
						Flag = "max  depth 不合格"
						#print("ref:" + refCiteBaseDic[Cite])
						#print(yesBaseLst)
						#print(yesBaseNumLst)
						#print(maxBase)
						#print(maxNum)
						#print(Alle)
						#print(AlleFreq)
						#print()
		
				else:
					if len(yesBaseLst) == 1:
						Alle = "NO"
						AlleFreq = "NO"
					else:
						del yesBaseLst[maxIndex]
						del yesBaseNumLst[maxIndex]
						#del All_StrandbiaDic[Sample][Cite].split("_")[maxIndex]

						secNum = max(yesBaseNumLst)
						secIndex = yesBaseNumLst.index(max(yesBaseNumLst))  # 求最大值的索引
						secBase = yesBaseLst[secIndex]	
						secStrand_bia = All_StrandbiaDic[Sample][Cite].split("_")[Strand_bia_index[secBase]]
						Strand_bia_Flag = secStrand_bia
						if TotNum >= DEP_THRES and secNum >= MIN_DEP_THRES:						
							if ( secStrand_bia != "na" and (float(secStrand_bia) >= STRANDED_RATIO_THRES ) and (float(secStrand_bia) <= 1-STRANDED_RATIO_THRES)  ):
								Alle = secBase
								AlleFreq = round(float(secNum)/TotNum,4)
								if AlleFreq < float(FREQ_THRES):
									AlleFreq == "NO"
									Alle = "NO"
							else:
								Alle = "NO"
								AlleFreq = "NO"
								Flag = "second  depth 合格  链偏性不合格"
							#	print("ref:" + refCiteBaseDic[Cite])
							#	print(yesBaseLst)
							#	print(yesBaseNumLst)
							#	print(maxBase)
							#	print(maxNum)
							#	print(Alle)
							#	print(AlleFreq)
							#	print()
	
						else:
							Alle = "NA"
							AlleFreq = "NA"
							Flag = "second  depth 不合格"
							#print("ref:" + refCiteBaseDic[Cite])
							#print(yesBaseLst)
							#print(yesBaseNumLst)
							#print(maxBase)
							#print(maxNum)
							#print(Alle)
							#print(AlleFreq)
							#print()
				AllAlle_Dic[Sample][Cite] = Alle
				AllAlleFreq_Dic[Sample][Cite] = AlleFreq

			else:
				Alle = "NA"
				AlleFreq = "NA"
			#if Cite == 45368:						
				#print("ref:" + refCiteBaseDic[Cite])
				#print(yesBaseLst)
				#print(yesBaseNumLst)
				#print(maxBase)
				#print(maxNum)
				#print(Alle)
				#print(AlleFreq)
				#print(Strand_bia)
				#print(Flag)
				#print(Strand_bia_Flag)
				#print()

	return AllAlle_Dic,AllAlleFreq_Dic









def outFreq(refCiteBaseDic,AllAlleFreq_Dic,FREQ_THRES,outputP):
	#outFreqlineDic = {}
	outFreqLst = []

	AllSampSNVCountDic = {}
	AllSampSNPCountDic = {}

	for SampID in AllAlleFreq_Dic:
		AllSampSNVCountDic[SampID] = 0
		AllSampSNPCountDic[SampID] = 0


	for refCite in refCiteBaseDic:
		Freq_line_Lst = [str(refCite)]
		NANO_count = 0
		SNVcount = 0
		SNPcount = 0
		freqLst = []
		for SampID in AllAlleFreq_Dic:
			if refCite in AllAlleFreq_Dic[SampID]:
				AlleFreq = AllAlleFreq_Dic[SampID][refCite]
			else:
				AlleFreq = "NA"
			freqLst.append(str(AlleFreq))
			if AlleFreq == "NA" or AlleFreq == "NO":
				NANO_count += 1
			else:
				if float(AlleFreq) >= FREQ_THRES and float(AlleFreq) < 1- FREQ_THRES:
					SNVcount += 1
					AllSampSNVCountDic[SampID] += 1
		

				elif AlleFreq >= FREQ_THRES:
					SNPcount += 1
					AllSampSNPCountDic[SampID] += 1
				else:
					NANO_count += 1


		Freq_line_Lst.append(str(SNVcount) + "/" + str(SNPcount))
		Freq_line_Lst.append("\t".join(freqLst))

		if NANO_count != len(AllAlleFreq_Dic):
			Freq_line = "\t".join(Freq_line_Lst)
			#outFreqlineDic[refCite] = Freq_line
			outFreqLst.append(Freq_line)


	out = "\n".join(outFreqLst)
	headLst = ["#Posi","snv/SNP"]
	for SampID in AllAlleFreq_Dic:
		headLst.append(SampID)
	header = "\t".join(headLst)

	outbigtable = outputP + "/all.iSNV_with_SNP.pyResults.txt"
	if ( os.path.exists(outbigtable)):
		os.remove(outbigtable)
	outAllsnvSNP_tableO = open(outbigtable,'a')
	outAllsnvSNP_tableO.write(header + "\n")
	outAllsnvSNP_tableO.write(out)
	outAllsnvSNP_tableO.close()



	return 	AllSampSNVCountDic,AllSampSNPCountDic





def outStat(Stat_Dic,AllSampSNVCountDic,AllSampSNPCountDic,outputP):
	outLst = []
	for sample in AllSampSNVCountDic:
		outlineLst = []
		outlineLst.append(sample)
		outlineLst.append(str(Stat_Dic["validPosNo"][sample]))
		outlineLst.append(str(Stat_Dic["IsSampleValid"][sample]))
		outlineLst.append(str(Stat_Dic["meanSiteDep"][sample]))
		outlineLst.append(str(AllSampSNVCountDic[sample]))
		outlineLst.append(str(AllSampSNPCountDic[sample]))
		outLst.append("\t".join(outlineLst))
	out = "\n".join(outLst)
	header = "\t".join(["Samp","validPosNo","IsSampleValid","meanSiteDep","no.SNV","np.SNP"])
	outStatF = outputP + "/samps.summary.txt"
	if ( os.path.exists(outStatF)):
		os.remove(outStatF)
	outStatFO = open(outStatF,'a')
	outStatFO.write(header + "\n")
	outStatFO.write(out)
	outStatFO.close()


def outStatofCov(refCiteBaseDic,All_CovDic):
	outLst = []
	for refCite in refCiteBaseDic:	
		outlineLst = [str(refCite)]
		NACount = 0
		for cov_sampID in All_CovDic:
			if refCite in All_CovDic[cov_sampID] :
				Cov = All_CovDic[cov_sampID][refCite]
			else:
				Cov = "NA"
				NACount += 1

			outlineLst.append(str(Cov))
		if NACount != len(All_CovDic):
			outLst.append("\t".join(outlineLst))
	out = "\n".join(outLst)
	HeaderLst = ["#Samp(A-G-C-T-tot)"]
	for cov_sampID in All_CovDic:
		HeaderLst.append(cov_sampID)
	header = "\t".join(HeaderLst)
	outStatF = outputP + "/out.samps-Coverage.bigtable.txt"
	if ( os.path.exists(outStatF)):
		os.remove(outStatF)
	outStatFO = open(outStatF,'a')
	outStatFO.write(header + "\n")
	outStatFO.write(out)
	outStatFO.close()




def outStatof_AGCT_Starnd_Bia(refCiteBaseDic,All_StrandbiaDic):
	outLst = []
	for refCite in refCiteBaseDic:	
		outlineLst = [str(refCite)]
		NACount = 0
		for cov_sampID in All_StrandbiaDic:
			if refCite in All_StrandbiaDic[cov_sampID] :
				Strand_bia_patt = All_StrandbiaDic[cov_sampID][refCite]
			else:
				Strand_bia_patt = "NA"
				NACount += 1 

			outlineLst.append(str(Strand_bia_patt))
		if NACount != len(All_StrandbiaDic):
			#print(NACount)
			#print(len(All_StrandbiaDic))
			outLst.append("\t".join(outlineLst))
	out = "\n".join(outLst)
	HeaderLst = ["#Samp"]
	for cov_sampID in All_StrandbiaDic:
		HeaderLst.append(cov_sampID)
	header = "\t".join(HeaderLst)
	outStatF = outputP + "/out.samps-AlleStrandBia.bigtable.txt"
	if ( os.path.exists(outStatF)):
		os.remove(outStatF)
	outStatFO = open(outStatF,'a')
	outStatFO.write(header + "\n")
	outStatFO.write(out)
	outStatFO.close()




def outStatof_ntpatt(refCiteBaseDic,All_ntpattDic):
	outLst = []
	for refCite in refCiteBaseDic:	
		outlineLst = [str(refCite)]
		NACount = 0
		for sampID in All_ntpattDic:
			#print(All_ntpattDic[sampID])
			if refCite in All_ntpattDic[sampID] :
				#print( All_ntpattDic[sampID])
				nt_patt = All_ntpattDic[sampID][refCite]
			else:
				nt_patt = "NA"
				NACount += 1
			outlineLst.append(str(nt_patt))
		if NACount != len(All_ntpattDic):
			outLst.append("\t".join(outlineLst))
	out = "\n".join(outLst)
	HeaderLst = ["#Samp"]
	for sampID in All_ntpattDic:
		HeaderLst.append(sampID)
	header = "\t".join(HeaderLst)
	outStatF = outputP + "/out.samps-ntfreq.bigtable.txt"
	if ( os.path.exists(outStatF)):
		os.remove(outStatF)
	outStatFO = open(outStatF,'a')
	outStatFO.write(header + "\n")
	outStatFO.write(out)
	outStatFO.close()


	



def main():
	ntfreqFileLst = readntfreqFList_configFile(ntfreqFList_configF)
	refID,refCiteBaseDic = refgenome_Base(refgenomeFile)

	print("total " + str(len(ntfreqFileLst)) + " ntfreq files")
	print()

	All_ntpattDic,Stat_Dic,All_StrandbiaDic,All_CovDic,All_ntfreqCitesLst = All_ntfreq_file(refID,ntfreqFileLst,DEP_THRES)

	##
	AllAlle_Dic,AllAlleFreq_Dic = find_alle(refCiteBaseDic,All_ntpattDic,All_StrandbiaDic,DEP_THRES,MIN_DEP_THRES,STRANDED_RATIO_THRES,FREQ_THRES)

	##output Freq bigtable
	AllSampSNVCountDic,AllSampSNPCountDic = outFreq(refCiteBaseDic,AllAlleFreq_Dic,FREQ_THRES,outputP)

	outStat(Stat_Dic,AllSampSNVCountDic,AllSampSNPCountDic,outputP)

	outStatofCov(refCiteBaseDic,All_CovDic)
	outStatof_AGCT_Starnd_Bia(refCiteBaseDic,All_StrandbiaDic)
	outStatof_ntpatt(refCiteBaseDic,All_ntpattDic)




if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)
	parser.add_option("-i","--ntfreqFList_configF",
	                  dest = "ntfreqFList_configF",
	                  default = "",
	                  metavar = "file",
	                  help = "Config File (all ntfreq File list,format of each line: path/sampID.ntfreq).  [required]")

	parser.add_option("-r","--refgenomeFile",
	                  dest = "refgenomeFile",
	                  default = "",
	                  metavar = "file",
	                  help = "Ref genome file.  [required]")

	parser.add_option("-o","--outputP",
	                  dest = "outputP",
	                  default = "",
	                  metavar = "path",
	                  help = "output path.  [required]")


	parser.add_option("-D","--DEP_THRES",
	                  dest = "DEP_THRES",
	                  default = "50",
	                  metavar = "int",
	                  help = "DEP_THRES.  [optional]")

	parser.add_option("-A","--MIN_DEP_THRES",
	                  dest = "MIN_DEP_THRES",
	                  default = "5",
	                  metavar = "int",
	                  help = "MIN_DEP_THRES.  [optional]")
	parser.add_option("-N","--VALID_SIZE_THRES",
	                  dest = "VALID_SIZE_THRES",
	                  default = "100",
	                  metavar = "int",
	                  help = "VALID_SIZE_THRES.  [optional]")
	parser.add_option("-F","--FREQ_THRES",
	                  dest = "FREQ_THRES",
	                  default = "0.02",
	                  metavar = "float",
	                  help = "FREQ_THRES.  [optional]")
	parser.add_option("-S","--STRANDED_RATIO_THRES",
	                  dest = "STRANDED_RATIO_THRES",
	                  default = "0.0",
	                  metavar = "float",
	                  help = "STRANDED_RATIO_THRES.  [optional]")

	(options,args) = parser.parse_args()
	ntfreqFList_configF     = os.path.abspath(options.ntfreqFList_configF)
	refgenomeFile     = os.path.abspath(options.refgenomeFile)
	outputP = os.path.abspath(options.outputP)
	DEP_THRES     = int(options.DEP_THRES)  ##最低深度
	MIN_DEP_THRES     = int(options.MIN_DEP_THRES)   ##支持alt的最小reads数
	VALID_SIZE_THRES     = int(options.VALID_SIZE_THRES)    ##样本有效需要的有效位点数
	FREQ_THRES     = float(options.FREQ_THRES)  ###alt的频率最小值
	STRANDED_RATIO_THRES     = float(options.STRANDED_RATIO_THRES)   ## reads链偏性检测



	main()




