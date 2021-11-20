
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
			seqID = InfoFl_Tags[18]
			clade = InfoFl_Tags[20]
			if seqID not in ["P1-S1","P1-S2","P1-S3","P2-S1","P5-S17","P7-S3"]:	
				if person in clade:
					if clade  not in cladeSampDic:
						cladeSampDic[clade] = []
					cladeSampDic[clade].append(seqID)
	return cladeSampDic





def tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,FiltPosiDic,repeatRegPosi):
	Sampsnv_Lst = {}
	SampSNP_Lst = {}
	FreqDic = samp_FreqDic[sample]
	tongji_Dic = {}

	snvCount,SNPCount = 0,0
	print(list(FreqDic.keys())[0:10])
	for Posi in FreqDic:
		if (Posi not in personNAPosi) and (Posi not in personCommonPosi) and (Posi not in FiltPosiDic)  and  (Posi not in repeatRegPosi):
			Freq = FreqDic[Posi]
			if Freq >= snv_MinFreq and Freq < snv_MaxFreq:
				snvCount += 1
				Sampsnv_Lst[Posi] = ''

			elif Freq >= snv_MaxFreq :
				SNPCount += 1
				SampSNP_Lst[Posi] = ''
	return Sampsnv_Lst,SampSNP_Lst




def readSnpEffF(sample,sampleSnpEffP,Sampsnv_Lst,SampSNP_Lst,CladeSNV_snpeffGeneDic,CladeSNP_snpeffGeneDic,snv_interNumDic,allEff_Stat_Dic,allFreqDic,PosiGeneDic):
	#PosiGeneDic = {}
	snpeffGeneDic = {}
	EffDic = {}
	CladeSNV_snpeffGeneDic[sample] = {}
	CladeSNP_snpeffGeneDic[sample] = {}
	allFreqDic[sample] = {}

	SNV_intraGeneNum = 0
	SNV_interGeneNum = 0
	SNP_intraGeneNum = 0
	SNP_interGeneNum = 0

	FreqDic = {}	
	sampleSnpEffF = sampleSnpEffP + "/" + sample + ".sort.rmdup.mpileup.minFreq0.05.snpeff.vcf"
	for sampleSnpEffFl in open(sampleSnpEffF).readlines():
		if "#" not in sampleSnpEffFl and sampleSnpEffFl != "\n":
			sampleSnpEffFl_Tags = sampleSnpEffFl.split("\n")[0].split("\t")
			Posi = int(sampleSnpEffFl_Tags[1])
			Infos = sampleSnpEffFl_Tags[7].split(";")
			
			if len(Infos) >5:
				Ann = Infos[5].split(",")[0].split("|")
				Gene = Ann[3]				
				Eff = Ann[1]
				EffDic[Posi] = Eff
				if Eff != "upstream_gene_variant" and Eff != "downstream_gene_variant":
					snpeffGeneDic[Posi] = Gene
					PosiGeneDic[Posi] = Gene
				else:
					if Posi in Sampsnv_Lst:
						SNV_interGeneNum += 1
					elif Posi in SampSNP_Lst:
						SNP_interGeneNum += 1
				#1/1:255:355:0:355:100.0%
				
				FreqInfo= sampleSnpEffFl_Tags[9]
				Freq = float(FreqInfo.split(":")[5].split("%")[0])
				

				FreqDic[Posi] = Freq



	snvSnpEffStatDic = {}
	SNPSnpEffStatDic = {}



	SNV_syn_Num,SNV_Nonsyn_Num,SNP_syn_Num,SNP_Nonsyn_Num = 0,0,0,0

	Eff_Stat_Dic = {}


	for Posi in snpeffGeneDic:
		GENE = snpeffGeneDic[Posi]
		Eff = EffDic[Posi]


		if Posi in Sampsnv_Lst:
			if GENE not in snvSnpEffStatDic:
				snvSnpEffStatDic[GENE] = 0
			snvSnpEffStatDic[GENE] += 1

		
			if Eff != "upstream_gene_variant" and Eff != "downstream_gene_variant":
				SNV_intraGeneNum += 1
				if Eff == "synonymous_variant":
					Flag = "SNV_syn"
					SNV_syn_Num += 1
					Eff_Flag = "iSNV_Syn"
				elif Eff != "synonymous_variant":
					Flag = "SNV_nonsyn"
					SNV_Nonsyn_Num += 1
					Eff_Flag = "iSNV_Miss"
					if Eff != "missense_variant":
						Eff_Flag = "iSNV_nonCoding"

				if GENE not in Eff_Stat_Dic:
					Eff_Stat_Dic[GENE] = {}
				if Flag  not in Eff_Stat_Dic[GENE]:
					Eff_Stat_Dic[GENE][Flag] = 0
				Eff_Stat_Dic[GENE][Flag] += 1

				allFreqDic[sample][Posi] = {}
				allFreqDic[sample][Posi][Eff_Flag] = FreqDic[Posi]




		elif Posi in SampSNP_Lst:
			if GENE not in SNPSnpEffStatDic:
				SNPSnpEffStatDic[GENE] = 0
			SNPSnpEffStatDic[GENE] += 1

			if Eff != "upstream_gene_variant" and Eff != "downstream_gene_variant":
				SNP_intraGeneNum += 1
				if Eff == "synonymous_variant":
					Flag = "SNP_syn"
					SNP_syn_Num += 1
				elif Eff != "synonymous_variant":
					Flag = "SNP_nonsyn"
					SNP_Nonsyn_Num += 1
				
				if GENE not in Eff_Stat_Dic:
					Eff_Stat_Dic[GENE] = {}
				if Flag  not in Eff_Stat_Dic[GENE]:
					Eff_Stat_Dic[GENE][Flag] = 0
				Eff_Stat_Dic[GENE][Flag] += 1




	#snv_interNumDic[sample] = [SNV_intraGeneNum,SNV_interGeneNum,SNP_intraGeneNum,SNP_interGeneNum,SNV_syn_Num,SNV_Nonsyn_Num,SNP_syn_Num,SNP_Nonsyn_Num]
	allEff_Stat_Dic[sample] = Eff_Stat_Dic
	CladeSNV_snpeffGeneDic[sample] = snvSnpEffStatDic
	CladeSNP_snpeffGeneDic[sample] = SNPSnpEffStatDic



	return CladeSNV_snpeffGeneDic,CladeSNP_snpeffGeneDic,snv_interNumDic,allEff_Stat_Dic,allFreqDic,PosiGeneDic



def FiltEcoliPosi(FiltEcoliPosiF):
	FiltPosiDic = {}
	for FiltEcoliPosiFl in open(FiltEcoliPosiF).readlines():
		if FiltEcoliPosiFl != "\n":
			FiltPosiDic[FiltEcoliPosiFl.split("\n")[0].split("\t")[0]] = ''

	return FiltPosiDic



def FiltRefRepeatRegions1(RepeatMaskedGnm):
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
	person = iSNV_SNP_table.split("/")[-2]
	if person == "P5_withP2":
		person = "P5"
	lieD,lsD = iSNV_SNP_tableRead(iSNV_SNP_table)

	FiltPosiDic = FiltEcoliPosi(FiltEcoliPosiF)
	repeatRegPosi = FiltRefRepeatRegions(RepeatMaskedGnm)

	print(len(repeatRegPosi))
	print(list(repeatRegPosi.keys())[0:10])
	print(list(lsD.keys())[0:10])
	print("hhhhh")
	cladeSampDic = readInfo(InfoF,person)
	sampleSnpEffP = snpEffPath
	#print(person)
	#print(cladeSampDic)
	FiltSamps = ["P1-S1","P1-S2","P1-S3","P2-S1","P5-S17","P7-S3"]
	

	for clade in cladeSampDic:
		#print(clade)
		cladeSamples = cladeSampDic[clade]
		
		cladeSamps = []
		for cladeSample in cladeSamples:
			for lieSamp in lieD:
				if lieSamp != "#Posi" and lieSamp != "snv/SNP" and cladeSample == lieSamp:
					cladeSamps.append(lieSamp)

		cladeSamps.sort()
		#print(cladeSamps)
		samp_FreqDic,sampNA_PosiDic = samp_Freq_Dic(cladeSamps,lieD,lsD)
		personNAPosi = Person_NA_Posi(cladeSamps,sampNA_PosiDic)
		personCommonPosi = Person_commonSNP_Posi(cladeSamps,samp_FreqDic,lsD)

		CladeSNV_snpeffGeneDic = {}
		CladeSNP_snpeffGeneDic = {}
		snv_interNumDic = {}
		allEff_Stat_Dic = {}
		allFreqDic = {}
		PosiGeneDic ={}

		for sample in cladeSamps:
			Sampsnv_Lst,SampSNP_Lst = tongjiFREQ(sample,samp_FreqDic,personNAPosi,personCommonPosi,FiltPosiDic,repeatRegPosi)		
			#print(Sampsnv_Lst["42319"])
			CladeSNV_snpeffGeneDic,CladeSNP_snpeffGeneDic,snv_interNumDic,allEff_Stat_Dic,allFreqDic,PosiGeneDic = readSnpEffF(sample,sampleSnpEffP,Sampsnv_Lst,SampSNP_Lst,CladeSNV_snpeffGeneDic,CladeSNP_snpeffGeneDic,snv_interNumDic,allEff_Stat_Dic,allFreqDic,PosiGeneDic)



		#print(PosiGeneDic)


##out SNV Eff Freq of each posi in each clade
		print("out SNV Eff Freq")
		Header = "\t".join(["Sample","posi","Flag","Freq"])
		outlines = []
		for sample in cladeSamps:
			for Posi in allFreqDic[sample]:
			 	for Eff_Flag in allFreqDic[sample][Posi]:
			 		if Eff_Flag != "iSNV_nonCoding":
			 	 		line = "\t".join([sample,str(Posi),Eff_Flag,str(allFreqDic[sample][Posi][Eff_Flag])])
			 	 		outlines.append(line)
		out = "\n".join(outlines)


		out_SNVFreq_F = out_P + "/" + clade + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".Eff_Freq.txt"
		if ( os.path.exists(out_SNVFreq_F)):
			os.remove(out_SNVFreq_F)
		out_SNP_F_O=open(out_SNVFreq_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(out + "\n")
		out_SNP_F_O.close()	



##out SNV Eff Freq of each gene in each clade
		GeneNumCount = {}
		for sample in cladeSamps:
			for Posi in allFreqDic[sample]:
			 	for Eff_Flag in allFreqDic[sample][Posi]:
			 		if Eff_Flag != "iSNV_nonCoding":
			 			geneID = PosiGeneDic[Posi]
			 			if geneID not in GeneNumCount :
			 				GeneNumCount[geneID]  = 0
			 			GeneNumCount[geneID] += 1
			 	 		


		print("out SNV Eff Freq")
		Header = "\t".join(["Gene","posi","Flag","Freq"])
		outlines = []
		for sample in cladeSamps:
			for Posi in allFreqDic[sample]:
			 	for Eff_Flag in allFreqDic[sample][Posi]:
			 		if Eff_Flag != "iSNV_nonCoding" and GeneNumCount[PosiGeneDic[Posi]] >= 3*len(cladeSamps) :
			 	 		line = "\t".join([PosiGeneDic[Posi],str(Posi),Eff_Flag,str(allFreqDic[sample][Posi][Eff_Flag])])
			 	 		outlines.append(line)
		out = "\n".join(outlines)


		out_SNVFreq_F = out_P + "/" + clade + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".gene.Eff_Freq.txt"
		if ( os.path.exists(out_SNVFreq_F)):
			os.remove(out_SNVFreq_F)
		out_SNP_F_O=open(out_SNVFreq_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(out + "\n")
		out_SNP_F_O.close()	





### out
		##person SNV snpeff Gene Stat
		PersonSNV_geneLst = []
		for sample in cladeSamps:
			for Gene in CladeSNV_snpeffGeneDic[sample]:
				#if CladeSNV_snpeffGeneDic[sample][Gene] > 1 and Gene not in PersonSNV_geneLst:
				if Gene not in PersonSNV_geneLst:
					PersonSNV_geneLst.append(Gene)

		linesLst = []
		for Gene in PersonSNV_geneLst:
			lineLst = [Gene]
			allNum = 0
			for sample in cladeSamps:
				if Gene in CladeSNV_snpeffGeneDic[sample]:
					Num = CladeSNV_snpeffGeneDic[sample][Gene]
				else:
					Num = 0
				allNum += Num

				lineLst.append(str(Num))
			lineLst.append(str(allNum))
			line = "\t".join(lineLst)
			if int(allNum) >= int(3 * len(cladeSamps)):
				linesLst.append(line)

		lines = "\n".join(linesLst)


		outHeader = ["gene"]
		for sample in cladeSamps:
			outHeader.append(sample)
		outHeader.append("SumNum")
		Header = "\t".join(outHeader)


		out_SNV_F = out_P + "/" + clade + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".snpEff_Stat.txt"
		if ( os.path.exists(out_SNV_F)):
			os.remove(out_SNV_F)
		out_SNP_F_O=open(out_SNV_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(lines + "\n")
		out_SNP_F_O.close()	

##-----------
		##person SNV  syn Stat
		Eff_PersonSNV_geneLst = []
		for sample in cladeSamps:
			for Gene in allEff_Stat_Dic[sample]:
				#if allEff_Stat_Dic[sample][Gene] > 1 and Gene not in PersonSNV_geneLst:
				if Gene not in Eff_PersonSNV_geneLst:
					Eff_PersonSNV_geneLst.append(Gene)


		linesLst = []
		for Gene in Eff_PersonSNV_geneLst:
			lineLst = [Gene]
			allNum = 0
			for sample in cladeSamps:
				if Gene in allEff_Stat_Dic[sample] and "SNV_syn" in allEff_Stat_Dic[sample][Gene]:
					Num = allEff_Stat_Dic[sample][Gene]["SNV_syn"]
				else:
					Num = 0
				allNum += Num

				lineLst.append(str(Num))
			lineLst.append(str(allNum))
			line = "\t".join(lineLst)
			if allNum >= 3 * len(cladeSamps):
				linesLst.append(line)


		lines = "\n".join(linesLst)



		outHeader = ["gene"]
		for sample in cladeSamps:
			outHeader.append(sample)
		outHeader.append("SumNum")
		Header = "\t".join(outHeader)


		out_SNV_F = out_P + "/" + clade + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".snpEff_Stat.Syn.txt"
		if ( os.path.exists(out_SNV_F)):
			os.remove(out_SNV_F)
		out_SNP_F_O=open(out_SNV_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(lines + "\n")
		out_SNP_F_O.close()	




######
		##person SNV  Nonsyn Stat

		Eff_PersonSNV_geneLst = []
		for sample in cladeSamps:
			for Gene in allEff_Stat_Dic[sample]:
				#if allEff_Stat_Dic[sample][Gene] > 1 and Gene not in PersonSNV_geneLst:
				if Gene not in Eff_PersonSNV_geneLst:
					Eff_PersonSNV_geneLst.append(Gene)


		linesLst = []
		for Gene in Eff_PersonSNV_geneLst:
			lineLst = [Gene]
			allNum = 0
			for sample in cladeSamps:
				if Gene in allEff_Stat_Dic[sample] and "SNV_nonsyn" in allEff_Stat_Dic[sample][Gene]:
					Num = allEff_Stat_Dic[sample][Gene]["SNV_nonsyn"]
				else:
					Num = 0
				allNum += Num

				lineLst.append(str(Num))
			lineLst.append(str(allNum))
			line = "\t".join(lineLst)
			if allNum >= 3 * len(cladeSamps):
				
				linesLst.append(line)


		lines = "\n".join(linesLst)



		outHeader = ["gene"]
		for sample in cladeSamps:
			outHeader.append(sample)
		outHeader.append("SumNum")
		Header = "\t".join(outHeader)


		out_SNV_F = out_P + "/" + clade + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".snpEff_Stat.NonSyn.txt"
		if ( os.path.exists(out_SNV_F)):
			os.remove(out_SNV_F)
		out_SNP_F_O=open(out_SNV_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(lines + "\n")
		out_SNP_F_O.close()	

		print()
		print()






'''
## snv Eff of each sample 
		Header = "\t".join(["Sample","posi","Flag","Freq"])	
		for sample in cladeSamps:
			outlines = []
			#print(sample)
			for Posi in allFreqDic[sample]:
			 	for Eff_Flag in allFreqDic[sample][Posi]:
			 	#for Eff_Flag in ["iSNV_Miss","iSNV_Syn"]:
			 		if Eff_Flag != "iSNV_nonCoding":
			 	 		line = "\t".join([sample,Posi,Eff_Flag,str(allFreqDic[sample][Posi][Eff_Flag])])
			 	 		outlines.append(line)
			out = "\n".join(outlines)

			#print(sample.split("_BDM")[0].split("BJ13-")[1]+ "hh") 

			out_SNVFreq_F = out_P + "/" + sample + ".snv_" + str(snv_MinFreq) + "-" +  str(snv_MaxFreq) + ".Eff_Freq.txt"
			if ( os.path.exists(out_SNVFreq_F)):
				os.remove(out_SNVFreq_F)
			out_SNP_F_O=open(out_SNVFreq_F,'a')
			out_SNP_F_O.write(Header + "\n")
			out_SNP_F_O.write(out + "\n")
			out_SNP_F_O.close()	
'''



'''
##person SNP snpeff Stat of SNP
		PersonSNP_geneLst = []
		for sample in cladeSamps:
			for Gene in CladeSNP_snpeffGeneDic[sample]:
				#if CladeSNP_snpeffGeneDic[sample][Gene] > 1 and Gene not in PersonSNP_geneLst:
				if Gene not in PersonSNP_geneLst:
					PersonSNP_geneLst.append(Gene)

		linesLst = []
		for Gene in PersonSNP_geneLst:
			lineLst = [Gene]
			allNum = 0
			for sample in cladeSamps:
				if Gene in CladeSNP_snpeffGeneDic[sample]:
					Num = CladeSNP_snpeffGeneDic[sample][Gene]
				else:
					Num = 0
				allNum += Num
				lineLst.append(str(Num))
			lineLst.append(str(allNum))
			line = "\t".join(lineLst)
			
			linesLst.append(line)

		lines = "\n".join(linesLst)

		out_SNP_F = out_P + "/" + clade + ".SNP_Over-" +  str(snv_MaxFreq) + ".snpEff_Stat.txt"
		if ( os.path.exists(out_SNP_F)):
			os.remove(out_SNP_F)
		out_SNP_F_O=open(out_SNP_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(lines + "\n")
		out_SNP_F_O.close()	







###---------------intra and inter gene
		#print(snv_interNumDic)

		PersonSNV_geneLst = []
		Header = "\t".join(["sample","SNV_intraGeneNum","SNV_interGeneNum","SNP_intraGeneNum","SNP_interGeneNum","SNV_syn_Num","SNV_Nonsyn_Num","SNP_syn_Num","SNP_Nonsyn_Num"])
		
		#SNV_intraNumLst = []
		#SNV_interNumLst = []
		outlines = []
		for sample in cladeSamps:
			#SNV_intraNumLst.append(snv_interNumDic[sample][0])
			#SNV_interNumLst.append(snv_interNumDic[sample][1])

			intra_inter_Line = [sample]
			for NUM in snv_interNumDic[sample]:
				intra_inter_Line.append(str(NUM))

			outlines.append("\t".join(intra_inter_Line))
		out = "\n".join(outlines)


		out_intra_inter_Gene_snvSNPStat_F = out_P + "/" + clade + ".snvSNP_" + str(snv_MinFreq) + "." +  str(snv_MaxFreq) + ".intra-inter.Gene.Stat.txt"
		if ( os.path.exists(out_intra_inter_Gene_snvSNPStat_F)):
			os.remove(out_intra_inter_Gene_snvSNPStat_F)
		out_SNP_F_O=open(out_intra_inter_Gene_snvSNPStat_F,'a')
		out_SNP_F_O.write(Header + "\n")
		out_SNP_F_O.write(out + "\n")
		out_SNP_F_O.close()	

'''




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

	parser.add_option("-E","--FiltEcoliPosiF",
					  dest = "FiltEcoliPosiF",
					  default = "",
					  metavar = "file",
					  help = "Filt Ecoli homo Posi file.  [required]")
	parser.add_option("-R","--RepeatMaskedGnm",
					  dest = "RepeatMaskedGnm",
					  default = "",
					  metavar = "file",
					  help = "RepeatMasker marked Gnm.  [required]")
	parser.add_option("-e","--snpEffPath",
					  dest = "snpEffPath",
					  default = "",
					  metavar = "path",
					  help = "snpEff file Path.  [required]")

	parser.add_option("-m","--snv_MinFreq",
					  dest = "snv_MinFreq",
					  default = "0.1",
					  metavar = "float",
					  help = "min Freq of snv to calculate (0-100%).  [required]")
	parser.add_option("-M","--snv_MaxFreq",
					  dest = "snv_MaxFreq",
					  default = "0.9",
					  metavar = "float",
					  help = "max Freq of snv to calculate (0-100%).  [required]")

	
	(options,args) = parser.parse_args()
	iSNV_SNP_table   = os.path.abspath(options.iSNV_SNP_table)
	out_P		  = os.path.abspath(options.out_P)
	InfoF   = os.path.abspath(options.InfoF)
	FiltEcoliPosiF   = os.path.abspath(options.FiltEcoliPosiF)
	RepeatMaskedGnm = os.path.abspath(options.RepeatMaskedGnm)
	snpEffPath = os.path.abspath(options.snpEffPath)
	
	snv_MinFreq		= float(options.snv_MinFreq)
	snv_MaxFreq		= float(options.snv_MaxFreq)

	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




