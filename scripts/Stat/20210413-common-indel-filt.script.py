
#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os
from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np 



def readInfo(InfoF):
	#cladeSampDic = {}
	allcladeSampDic = {}
	for InfoFl in open(InfoF).readlines():
		if "#" not in InfoFl:
			InfoFl_Tags = InfoFl.split("\n")[0].split("\t")
			seqID = InfoFl_Tags[1]
			clade = InfoFl_Tags[2]
			if clade  not in allcladeSampDic:
				allcladeSampDic[clade] = []
			allcladeSampDic[clade].append(seqID)

	return allcladeSampDic




def readPindelSV(pindelDelFs,allSV_homoDic,allnotHomoDic,AllHomoLen):
	pindelDic = {}	
	for samppindelSVF in pindelDelFs:
		print("del file: " + samppindelSVF)
		ID = samppindelSVF.split(".")[0]
		file = pindelDelP + "/" + samppindelSVF
		pindelDic[ID] = {"DEL":{}}
		allSV_homoDic[ID],allnotHomoDic[ID],AllHomoLen[ID] = {},{},{}
		for samppindelSVFl in open(file).readlines():
			if "#" not in samppindelSVFl and samppindelSVFl != "\n":
				Tags = samppindelSVFl.split("\t")
				infoTags = Tags[7].split(";")
				Posi = int(Tags[1])
				SVType = infoTags[3]
				
				refDp = int(Tags[9].split(":")[1].split(",")[0])
				altDp = int(Tags[9].split(":")[1].split(",")[1])
				Freq = round(float(altDp)/(refDp + altDp),2)
				for infTag in infoTags:
					if "SVTYPE" in infTag:
						SVType = infTag.split("=")[-1]

					if "SVLEN" in infTag:
						SVLen = abs(int(infTag.split("=")[-1]))
					if "HOMLEN" in infTag:
						HOMLEN = int(infTag.split("=")[1])
					if "HOMSEQ" in infTag:
						HOMSEQ = infTag.split("=")[1]
				
				if SVType == "DEL" and SVLen <= maxIndelLen and int(altDp) >= minAltNum and Freq >= minFreq:
					ref = Tags[3]
					alt = Tags[4]
					SV = ''
					for ii in range(len(alt),len(ref)):
						SV += ref[ii]
					pindelDic[ID][SVType][Posi] = SV + "_" + str(Freq)
					
					if HOMLEN >= 2:
						hmbs = {}
						for bs in HOMSEQ:
							if bs not in hmbs:
								hmbs[bs] = 0
							hmbs[bs] += 1
						if len(hmbs) == 1:
							allSV_homoDic[ID][Posi] = HOMSEQ
							AllHomoLen[ID][Posi] = HOMLEN
						else:
							allnotHomoDic[ID][Posi] = ''
					else:
						allnotHomoDic[ID][Posi] = ''

				
						
	return pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen




def readPindelIns(pindelInstFs,pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen):
	print(pindelDic.keys())
	for pindelInstF in pindelInstFs:
		print("ins file: " + pindelInstF)
		ID = pindelInstF.split(".")[0]
		file = pindelInsP + "/" + pindelInstF
		pindelDic[ID]["INS"] = {}
		#allSV_homoDic[ID] = {}
		for samppindelSVFl in open(file).readlines():
			if "#" not in samppindelSVFl and samppindelSVFl != "\n":
				Tags = samppindelSVFl.split("\t")
				infoTags = Tags[7].split(";")
				Posi = int(Tags[1])
				SVType = infoTags[3]
				refDp = int(Tags[9].split(":")[1].split(",")[0])
				altDp = int(Tags[9].split(":")[1].split(",")[1])
				Freq = round(float(altDp)/(refDp + altDp),2)

				for infTag in infoTags:
					if "SVTYPE" in infTag:
						SVType = infTag.split("=")[-1]
					if "SVLEN" in infTag:
						SVLen = abs(int(infTag.split("=")[-1]))
					if "HOMLEN" in infTag:
						HOMLEN = int(infTag.split("=")[1])
					if "HOMSEQ" in infTag:
						HOMSEQ = infTag.split("=")[1]
			
				if SVType == "INS" and SVLen <= maxIndelLen and int(altDp) >= minAltNum and Freq >= minFreq:
					ref = Tags[3]
					alt = Tags[4]
					SV = ''
					for ii in range(len(ref),len(alt)):
						SV += alt[ii]
					pindelDic[ID][SVType][Posi] = SV + "_" + str(Freq)
					if HOMLEN >= 2:
						hmbs = {}
						for bs in HOMSEQ:
							if bs not in hmbs:
								hmbs[bs] = 0
							hmbs[bs] += 1
						if len(hmbs) == 1:
							allSV_homoDic[ID][Posi] = HOMSEQ
							AllHomoLen[ID][Posi] = HOMLEN
						else:
							allnotHomoDic[ID][Posi] = ''

					else:
						allnotHomoDic[ID][Posi] = ''


	return pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen






def readVarscanSV(varscanFs):
	allSamp_VarscanDic = {}
	#print(varscanFs)
	for varscanF in varscanFs:
		#print(varscanF)
		ID = varscanF.split("/")[-1].split(".")[0]
		varscanF = varscanP + "/"+ varscanF	
		allSamp_VarscanDic[ID] = {"DEL":{},"INS":{}}
		for varscanFl in open(varscanF).readlines():
			if "Chrom" not in varscanFl and varscanFl != "\n":
				Tags = varscanFl.split("\t")
				
				Posi = int(Tags[1])
				if Posi != 1 :
					PASS = Tags[6]
					Freq = float(Tags[6].split("%")[0])*0.01
					
					read2Num = int(Tags[5])
					ref = Tags[2] 
					svtype = Tags[3].split("/")[1]
					if Freq >= minFreq and read2Num >= minAltNum:
						if "-" in svtype:
							SVType = "DEL"
							SV = svtype.split("-")[1]
						if "+" in svtype:
							SVType = "INS"
							SV = svtype.split("+")[1]
						if len(SV) <= maxIndelLen:
							allSamp_VarscanDic[ID][SVType][Posi] = SV + "_" + str(Freq)

	return allSamp_VarscanDic




def read_file_name(file_dir):
	Files = []
	for root,dirs,files in os.walk(file_dir):
		for file in files:
			if "_D.vcf" in file:
				Files.append(os.path.join(root,file))
	return Files






def sameSVofpindelAndVarscan(CladeSamps,pindelDic,allSamp_VarscanDic,sampReliabSVs,allnotCommom):
	AllcommSVPosiLst = []
	cladeSVDic = {}
	#print(allSamp_VarscanDic['P1-S14']["INS"])
	#print(pindelDic['P1-S14']["INS"])
	for sampleID in CladeSamps:
		sampReliabSVs[sampleID] = {}
		for SVType in pindelDic[sampleID]:	
			if SVType not in cladeSVDic:
				cladeSVDic[SVType] = {}		
			sampReliabSVs[sampleID][SVType] = {}
			for POS in pindelDic[sampleID][SVType]:
				pindelSV_freq = float(pindelDic[sampleID][SVType][POS].split("_")[1])
				SV = pindelDic[sampleID][SVType][POS].split("_")[0]
				varscanSV_freq = 0

				if POS in allSamp_VarscanDic[sampleID][SVType]:
					varscanSV_freq = float(allSamp_VarscanDic[sampleID][SVType][POS].split("_")[1])
					
				if varscanSV_freq != 0 and (pindelSV_freq  >= 0.9 or varscanSV_freq > 0.9) :
					sampReliabSVs[sampleID][SVType][POS] = SV 
					if POS == 5446395:
						print(sampleID)
						print(pindelSV_freq)
						print(varscanSV_freq)
						print()
					if POS not in cladeSVDic[SVType]:
						cladeSVDic[SVType][POS] = 0
					cladeSVDic[SVType][POS] += 1
				else:	
					if POS == 5446395:
						print(sampleID)

				#elif varscanSV_freq != 0 and (varscanSV_freq > 0.8) :
					#sampReliabSVs[sampleID][SVType][POS] = allSamp_VarscanDic[sampleID][SVType][POS].split("_")[0]
					#if POS not in cladeSVDic[SVType]:
						#cladeSVDic[SVType][POS] = 0
					#cladeSVDic[SVType][POS] += 1


	#print(cladeSVDic["INS"][5446395])
	#print("lalalal")
	for Type in cladeSVDic: 
		allnotCommom[Type]=[]
		for POSI in cladeSVDic[Type]:
			shareCount = cladeSVDic[Type][POSI]

			if shareCount != 0 and shareCount != len(CladeSamps):
				if POSI not in allnotCommom[Type]:
					allnotCommom[Type].append(POSI)




	return allnotCommom,sampReliabSVs

			
			




def main():

	pindelDelFs = os.listdir(pindelDelP)
	pindelInstFs = os.listdir(pindelInsP)
	varscanFs = os.listdir(varscanP)

	allSV_homoDic,allnotHomoDic,AllHomoLen = {},{},{}
	#print(allSV_homoDic)

	pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen = readPindelSV(pindelDelFs,allSV_homoDic,allnotHomoDic,AllHomoLen)
	pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen = readPindelIns(pindelInstFs,pindelDic,allSV_homoDic,allnotHomoDic,AllHomoLen)
	allSamp_VarscanDic = readVarscanSV(varscanFs)



	allcladeSampDic = readInfo(InfoF)


	allSamps = []	
	for Samp_delly in pindelDelFs:
		sampID = Samp_delly.split("/")[-1]
		allSamps.append(sampID)



	for clade in allcladeSampDic:
	#for clade in ["P4"]:
		#print(clade)
		cladeSamps = allcladeSampDic[clade]	
		hpySamps = []
		if clade == "P4" or clade == "P3_1":
			for samp in cladeSamps:
				if samp not in ["P4-S1","P4-S2","P4-S4","P4-S5"]:
					hpySamps.append(samp)

		P4Samps = []
		if clade == "P4":
			for samp in cladeSamps:
				if samp not in ["P4-S1","P4-S2","P4-S4","P4-S5"]:
					P4Samps.append(samp)
		P3_1_Samps = []
		if clade == "P3_1":
			for samp in cladeSamps:
				P3_1_Samps.append(samp)


		sampReliabSVs,allnotCommom = {},{}
		allnotCommom,sampReliabSVs = sameSVofpindelAndVarscan(cladeSamps,pindelDic,allSamp_VarscanDic,sampReliabSVs,allnotCommom)

		print(allnotCommom)
		print(len(allnotCommom["DEL"]))
		print(len(allnotCommom["INS"]))
		outDels = {}
		outDels_Pos = {}
		sampsHomo = []
		for samp in cladeSamps:
			count,outDels[samp] = {},{}
			for Type in ["DEL","INS"]:
				count[Type],outDels[samp][Type],outDels_Pos[Type] = {},{},[]
				for POSI in allnotCommom[Type]:
					SV = ''
					flag = ''
					if POSI in allSV_homoDic[samp]:
						if POSI in sampReliabSVs[samp][Type]:
							SV = sampReliabSVs[samp][Type][POSI]
							flag="Homo"
						else:
							flag="NotHas"
					elif POSI in allnotHomoDic[samp]:
						if POSI in sampReliabSVs[samp][Type]:
							flag="NotHomo"
							SV = "NotHomo"
						#if POSI in sampReliabSVs[samp][Type]:
							#SV = sampReliabSVs[samp][Type][POSI]
					#else:
						#print(sampReliabSVs[samp][Type])
						#print(samp)
						#print(POSI)
						#print()

					outDels[samp][Type][POSI] = SV
					if POSI not in outDels_Pos[Type]:
						outDels_Pos[Type].append(POSI)

				
					if flag not in count[Type]:
						count[Type][flag] = 0
					count[Type][flag] += 1
			
			if samp in hpySamps:
				hpyFlag = "hpy"
			else:
				hpyFlag = "nomor"
			
			if samp in P4Samps:
				hpyClade = "P4-hpy"
			elif samp in P3_1_Samps:
				hpyClade = "P3-1"
			else:
				hpyClade = "nomor"
		
			head = ["samp","hpyFlag","hpyClade"]
			l = [samp,hpyFlag,hpyClade]
			sumNum,Homosum,NotHomosum =0,0,0
			
			for flg in ["Homo","NotHomo"]:  #,"NotHas"
				for tp in ["DEL","INS"]:
					head.append(tp + "_" + flg)
					num = 0
					if tp in count and flg in count[tp]:
						num = count[tp][flg]
					sumNum += num
					if flg == "Homo":
						Homosum += num
					if flg == "NotHomo":
						NotHomosum += num

					l.append(str(num))	
			l.append(str(Homosum))
			l.append(str(NotHomosum))
			if sumNum != 0:
				per_Hom = str(round(float(Homosum)*100/sumNum,2))
				per_noHom = str(round(float(NotHomosum)*100/sumNum,2))
			else:
				per_Hom = "0"
				per_noHom = "0"
			l.append(per_Hom)
			l.append(per_noHom)
					
			l.append(str(sumNum))
			head.append("HomosumNum")
			head.append("NonHomsumNum")
			head.append("perct_Homo")
			head.append("perct_noHomo")

			head.append("sumNum")
			H = "\t".join(head)
			sampsHomo.append("\t".join(l))
		print(clade)
		print(H)
		print("\n".join(sampsHomo))
		print()
		out_F = out_P  + "/" +  clade + ".reliabIndels.Homo.notHomo.stat.txt"
		if ( os.path.exists(out_F)):
			os.remove(out_F)
		out_SNP_F_O=open(out_F,'a')
		out_SNP_F_O.write(H + "\n")
		out_SNP_F_O.write("\n".join(sampsHomo) + "\n")
		out_SNP_F_O.close()	


		out_F = out_P  + "/allClades.reliabIndels.Homo.notHomo.stat.txt"
		if ( not os.path.exists(out_F)):
			EFlg = "not"
		else:
			EFlg = ''
		
		out_SNP_F_O=open(out_F,'a')
		if EFlg == "not":
			out_SNP_F_O.write(H + "\n")
		out_SNP_F_O.write("\n".join(sampsHomo) + "\n")
		out_SNP_F_O.close()	


		#print(outDels["P1-S1"]["DEL"])
		H = ["typy","Posi","homoLen"]
		for samp in cladeSamps:
			H.append(samp)
		hed = "\t".join(H)
		for Type in ["DEL","INS"]:
			Ls = []
			for POSI in allnotCommom[Type]:
				for samp in cladeSamps:
					if samp in AllHomoLen and POSI in AllHomoLen[samp]:
						Homo_len = AllHomoLen[samp][POSI]
						break
				L=[Type,str(POSI),str(Homo_len)]

				for samp in cladeSamps:
					L.append(outDels[samp][Type][POSI])

				Ls.append("\t".join(L))
			out = "\n".join(Ls)
			#print(hed)
			#print(out)
			#print()

			out_F = out_P  + "/" +  clade  + "." + Type +  ".reliabIndels.Homo.notHomo.Posi.txt"
			if ( os.path.exists(out_F)):
				os.remove(out_F)
			out_SNP_F_O=open(out_F,'a')
			out_SNP_F_O.write(hed + "\n")
			out_SNP_F_O.write(out + "\n")
			out_SNP_F_O.close()	



'''
	print(sampReliabSVs)
	header = ["Posi"]
	allSamps.sort()
	for samp in allSamps:
		header.append(samp)
	Head = "\t".join(header)	

	outL = []
	#allSV_over50notCommom.sort()
	print(allSV_over50notCommom)
	for SV_over50notCommom in allSV_over50notCommom:
		outl = [str(SV_over50notCommom)]
		for samp in allSamps:
			Len = ''
			if SV_over50notCommom in sampReliabSVs[samp]:
				Len = str(sampReliabSVs[samp][SV_over50notCommom])
			outl.append(Len)
		print(outl)
		outL.append("\t".join(outl))

	out_F = out_P + "/allPerson.Deletion.pindel_delly_common.Stat.txt"
	if ( os.path.exists(out_F)):
		os.remove(out_F)
	out_SNP_F_O=open(out_F,'a')
	out_SNP_F_O.write(Head + "\n")
	out_SNP_F_O.write("\n".join(outL) + "\n")
	out_SNP_F_O.close()	

'''



if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-f","--InfoF",
					  dest = "InfoF",
					  default = "",
					  metavar = "file",
					  help = "Info File.  [required]")

	parser.add_option("-d","--pindelDelP",
					  dest = "pindelDelP",
					  default = "",
					  metavar = "path",
					  help = "pindel Del Path.  [required]")
	parser.add_option("-i","--pindelInsP",
					  dest = "pindelInsP",
					  default = "",
					  metavar = "path",
					  help = "pindel INS Path.  [required]")

	parser.add_option("-v","--varscanP",
					  dest = "varscanP",
					  default = "",
					  metavar = "path",
					  help = "varscanP.  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "out_P. (.common-) [required]")

	parser.add_option("-m","--maxIndelLen",
					  dest = "maxIndelLen",
					  default = "50",
					  metavar = "int",
					  help = "max indels length.  [optional]")

	parser.add_option("-a","--minFreq",
					  dest = "minFreq",
					  default = "0.7",
					  metavar = "int",
					  help = " min freq of alt.  [optional]")

	parser.add_option("-b","--minAltNum",
					  dest = "minAltNum",
					  default = "10",
					  metavar = "int",
					  help = "min alt num .  [optional]")



	(options,args) = parser.parse_args()
	InfoF   = os.path.abspath(options.InfoF)
	pindelDelP = os.path.abspath(options.pindelDelP)
	pindelInsP = os.path.abspath(options.pindelInsP)
	varscanP  = os.path.abspath(options.varscanP)
	out_P  = os.path.abspath(options.out_P)
	maxIndelLen =int(options.maxIndelLen)
	minFreq =float(options.minFreq)
	minAltNum =int(options.minAltNum)



	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)

	out_F = out_P  + "/allClades.reliabIndels.Homo.notHomo.stat.txt"
	if ( os.path.exists(out_F)):
		os.remove(out_F)

	main()





