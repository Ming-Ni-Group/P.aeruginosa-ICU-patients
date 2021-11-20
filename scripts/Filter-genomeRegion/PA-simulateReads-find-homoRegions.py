#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os
import sys
from optparse import OptionParser

import pysam



def readbam(bamF):
	linecount = 0
	PosiDic = {}
	bf = pysam.AlignmentFile(bamF, 'rb')
	for bam_line in bf:
		POS = bam_line.pos 
		CIGAR = bam_line.cigar
		#SEQ = bam_line.seq
		QNAME = bam_line.qname
		#SEQ_len = len(SEQ)	
		MAPQ =bam_line.mapq
		
		if MAPQ >= 0:	
			PAstart=int(QNAME.split("_")[1].split("-")[0])
			for CIGAR_grp in CIGAR:
				if CIGAR_grp[0] != 0:
					PAstart += CIGAR_grp[1]
				if CIGAR_grp[0] == 0 and CIGAR_grp[1] >= 30:
					matchCount = CIGAR_grp[1]
					
					break
			

			for indx in range(0,matchCount):
				if str(PAstart + indx) not in PosiDic:
					PosiDic[str(PAstart + indx)] = ''

		linecount += 1
		#if linecount > 100:
			#break
		if linecount%1000 == 0:
			print(linecount)
	return PosiDic






def main():
	PosiLst = readbam(bamF)
	print(len(PosiLst))
	lineL = []
	for Posi in PosiLst:
		line = "\t".join([Posi,refID])
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

	parser.add_option("-i","--bamF",
					  dest = "bamF",
					  default = "",
					  metavar = "path",
					  help = "bam file of .  [required]")
	parser.add_option("-o","--outputfile",
					  dest = "outputfile",
					  default = "",
					  metavar = "file",
					  help = "output path [required]")

	parser.add_option("-r","--refID",
					  dest = "refID",
					  default = "",
					  metavar = "str",
					  help = "refID [required]")



	(options,args) = parser.parse_args()
	bamF	 = os.path.abspath(options.bamF)
	outputfile	 = os.path.abspath(options.outputfile)
	refID	 = options.refID
	


	if (os.path.exists(outputfile)):
		os.remove(outputfile)

	

	main()

