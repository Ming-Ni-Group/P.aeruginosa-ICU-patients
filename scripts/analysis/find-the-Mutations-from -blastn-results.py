##
##This script is to search mutations that occur in multiple MLST isolates in this study 

#Help for this script
# python  screen-geneMutation.py -h 

## parameters 
# -g genes file of the reference 
# -m the recurrently mutations of the genes 
# -b blast results 
# -o output path of the file

## example 
# python find-the-Mutations-from -blastn-results.py -g ../../meta_data/PAO1.sequence-genes.txt  -m  ../../meta_data/list of recurrent genes and mutations recurrently.txt
# -o   $outpath  -b ../../meta_data


#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,re
from optparse import OptionParser


def readGene(GeneF):
	IDseqs = {}
	for GeneFl in open(GeneF).readlines():
		if ">" in GeneFl:
			ID = GeneFl.split("\n")[0].split(" ")
			idd = ID[1].split("=")[1].split("]")[0]
			if idd not in IDseqs:
				IDseqs[idd] = ''
		else:
			IDseqs[idd] += GeneFl.split("\n")[0]
	
	return IDseqs



def readMutation(MutationsF):
	IDs = {}
	for GeneFl in open(MutationsF).readlines():
		gene = GeneFl.split("\n")[0].split("\t")[0]
		geneMutat = GeneFl.split("\n")[0].split("\t")[1]
		posi = re.findall(r"\d+",geneMutat)[0]
		refBs= geneMutat.split(posi)[1].split(">")[0]
		alleBs= geneMutat.split(posi)[1].split(">")[1]
		print(posi)
		if "-" not in geneMutat:
			if gene not in IDs:
				IDs[gene] = {}
			IDs[gene][int(posi)] = refBs + "-" + alleBs 

	return IDs

def readblastResult(gene,blastResultP,mutats):
	print(gene)
	blastResultF = blastResultP + "/" + gene + "-Alignment.txt"
	chrom_query_bases = {}
	chrom_Sbjct_bases = {}
	chromStart = {}
	chrom = ''
	for blastResultFl in open(blastResultF).readlines():		
		#print(blastResultFl)
		if ">" in blastResultFl:
			chrom = blastResultFl.split("\n")[0].split(">")[1]
			if chrom not in chrom_query_bases:
				chrom_query_bases[chrom] = ''
			if chrom not in chrom_Sbjct_bases:
				chrom_Sbjct_bases[chrom] = ''
			if chrom not in chromStart:
				chromStart[chrom] = 0

		#print(chrom)
		if chrom != '':
			if "Query" in blastResultFl:
				#print(blastResultFl)

				if chrom_query_bases[chrom] == "":
					Start = re.findall(r"\d+",blastResultFl)[0]
					chromStart[chrom] = int(Start)
				baseidxs = []
				idx = 0
				for base in blastResultFl:
					if base in ["A","G","C","T","-","a","g","c","t"]:
						baseidxs.append(idx)
					idx += 1

	
				query = blastResultFl[baseidxs[0]:baseidxs[-1]+1 ]
				chrom_query_bases[chrom] += query
			if "Sbjct" in blastResultFl:
				
				#print(blastResultFl)
				#print("Sbjct")
				sbjct = blastResultFl[baseidxs[0]:baseidxs[-1]+1 ]
				chrom_Sbjct_bases[chrom] += sbjct


	#print(chrom_query_bases["Pseudomonas aeruginosa strain SE5443 chromosome, complete genome"])
	#print(chrom_Sbjct_bases["Pseudomonas aeruginosa strain SE5443 chromosome, complete genome"])

	outs = []
	for chromID in chrom_query_bases:
		blastseq_query = chrom_query_bases[chromID]
		blastseq_Sbjct = chrom_Sbjct_bases[chromID]
		cout = chromStart[chromID]-1

		idx = 0
		muts = mutats[gene]

		#print(blastseq_query[787])
		#print(blastseq_Sbjct)
		for bs in blastseq_query:
			
			if bs != "-":
				cout += 1
				if cout in muts:
					#print(chromID)
					#print(cout)
					#print(muts)
					ref_alle = muts[cout]
					ref = ref_alle.split("-")[0]
					alle = ref_alle.split("-")[1]
					#print(ref)
					#print(ref)
					#rawRef = blastseq_query[idx]
					sbjctBase = blastseq_Sbjct[idx]


					if alle == sbjctBase:
						
						l = []
						l.append(gene)
						#l.append(chromID.split(" ")[3])
						l.append(chromID)
						l.append(str(cout))
						l.append(ref_alle)
						outs.append("\t".join(l))
			idx += 1
	Os = "\n".join(outs)

	return Os




	


def main():
	IDseqs = readGene(GeneF)
	IDs = readMutation(MutationsF)
	outF = out_P + "/mutation.blastn.txt"
	if ( os.path.exists(outF)):
		os.remove(outF)

	for idd in IDs:		
		genSeq = IDseqs[idd]
		#if idd == "fptA":
		Os = readblastResult(idd,blastResultP,IDs)
		print(Os)	
		out_SNP_F_O=open(outF,'a')
		out_SNP_F_O.write(Os + "\n")
		out_SNP_F_O.close()		






if __name__ == "__main__":
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-g","--GeneF",
					  dest = "GeneF",
					  default = "",
					  metavar = "file",
					  help = "ref genome file .  [required]")

	parser.add_option("-o","--out_P",
					  dest = "out_P",
					  default = "",
					  metavar = "path",
					  help = "output stat file Path of snv Freq distribution.  [required]")

	parser.add_option("-m","--MutationsF",
					  dest = "MutationsF",
					  default = "",
					  metavar = "path",
					  help = "Region File.  [required]")

	parser.add_option("-b","--blastResultP",
					  dest = "blastResultP",
					  default = "",
					  metavar = "path",
					  help = "blast Result P.  [required]")




	(options,args) = parser.parse_args()
	GeneF   = os.path.abspath(options.GeneF)
	out_P		 = os.path.abspath(options.out_P)
	MutationsF   = os.path.abspath(options.MutationsF)
	blastResultP = os.path.abspath(options.blastResultP)
	if ( not os.path.exists(out_P)):
		os.mkdir(out_P)


	main()




