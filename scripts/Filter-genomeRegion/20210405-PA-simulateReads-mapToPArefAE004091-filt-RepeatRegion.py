#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os
import sys
from optparse import OptionParser

import pysam



def readBlastF1(blastF):
	PosiDic = {}
	for blastFl in open(blastF).readlines():
		Tags = blastFl.split("\n")[0].split("\t")
		seqid = Tags[0]
		indentity = float(Tags[2])
		qIdS = int(seqid.split("_")[1].split("-")[0])
		qS_index = int(Tags[6])
		qE_index = int(Tags[7])
		rS = int(Tags[8])
		rE = int(Tags[9])

		q_S = min([qS_index,qE_index])
		q_E = max([qS_index,qE_index])
		r_S = min([rS,rE])
		r_E = max([rS,rE])
		
		if seqid != "Posi_" + str(rS) + "-" + str(rE) and (q_E - q_S >= 30) and indentity >= 60 :
			for rPos in range(r_S,r_E+1):
				if rPos not in PosiDic:
					PosiDic[rPos] = ''
			for qPindex in range(q_S,q_E +1):
				if qIdS + qPindex -1 not in PosiDic:
					PosiDic[qIdS + qPindex -1] = ''

	return PosiDic





def readBlastF(blastF):
	PosiDic = {}
	seqqBlastCount = {}
	for blastFl in open(blastF).readlines():
		Tags = blastFl.split("\n")[0].split("\t")
		seqid = Tags[0]
		indentity = float(Tags[2])
		qIdS = int(seqid.split("_")[1].split("-")[0])
		qIdE = int(seqid.split("_")[1].split("-")[1])
		qS_index = int(Tags[6])
		qE_index = int(Tags[7])
		rS = int(Tags[8])
		rE = int(Tags[9])

		q_S = min([qS_index,qE_index])
		q_E = max([qS_index,qE_index])
		r_S = min([rS,rE])
		r_E = max([rS,rE])

		if seqid != "Posi_" + str(rS) + "-" + str(rE) and (q_E - q_S >= 30) and indentity >= 60 :
			for rPos in range(r_S,r_E+1):
				if rPos not in PosiDic:
					PosiDic[rPos] = ''
			for qPindex in range(qIdS,qIdE +1):
				PosiDic[qPindex] = ''

	return PosiDic








def main():
	PosiDic = readBlastF(blastF)
	PosiLst = []
	print(len(PosiDic))
	for Po in PosiDic:
		PosiLst.append(Po)
	PosiLst.sort()
	lineL = []
	for Posi in PosiLst:
		line = "\t".join([str(Posi),refID])
		lineL.append(line)
	lines = "\n".join(lineL)
	outFO = open(outputfile,'a')
	outFO.write(lines + "\n")
	outFO.close()




if __name__ == "__main__":
################################################################################################################################
# Parameters
################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--blastF",
					  dest = "blastF",
					  default = "",
					  metavar = "file",
					  help = "bam file of .  [required]")
	parser.add_option("-o","--outputfile",
					  dest = "outputfile",
					  default = "",
					  metavar = "file",
					  help = "output file [required]")

	parser.add_option("-r","--refID",
					  dest = "refID",
					  default = "AE004091",
					  metavar = "str",
					  help = "refID [required]")



	(options,args) = parser.parse_args()
	blastF	 = os.path.abspath(options.blastF)
	outputfile	 = os.path.abspath(options.outputfile)
	refID	 = options.refID
	


	if (os.path.exists(outputfile)):
		os.remove(outputfile)

	

	main()

