# Whole-genome sequence characterization of longitudinal Pseudomonas aeruginosa isolates from long-term inpatients in intense care units


## summary

Pseudomonas aeruginosa (P. aeruginosa) is a common opportunistic pathogen that is one of the leading causes of nosocomial infection. Many previous studies regarding the surveillance of nosocomial transmission and epidemiology of P. aeruginosa have been reported. However, longitudinal studies regarding the dynamics of P. aeruginosa infection in healthcare settings were rare. We obtained longitudinal samples from 7 elderly patients (P1 to P7) long-term cared in intensive care units (ICUs) in periods of 4 – 19 months. P. aeruginosa was isolated from all the patients, and whole-genome sequence analyses indicated ten phylogenetic distant clones (clone complex) from the 7 patients. In the sequential sputum samples of two patients (P3 and P5), turnover of the main clones were observed. One clone was identified in two patients, it appeared as a transient clone that frequent emerged and lost. We compared the intraclonal genomic diversity of these P. aeruginosa complex, and found hypermutator isolates from two clones, which had remarkably more single nucleotide polymorphisms and variations in homopolymers sequences than the others. Hypermutator clones were associated with mutations T47I-G521S and P27L in MutL protein. We also found 24 recurrent mutated genes that had intraclonal diversity in >= 2 clones, 8 of which involve in siderophore biosynthesis, antibiotic resistance and biofilm formation. Especially, one recurrent mutation, S698F in FtpA, was observed in four clone. Our study emphasizes the importance of longitudinal surveillance of nosocomial P. aeruginosa clones in the healthcare settings, and suggests convergent microevolution and adaption of P. aeruginosa infection within long-term ICU patients.




## Dependencies

- python 3.7
- perl version v5.26.2
- RAxML 8.2.12
- R version 3.5.1
- BWA  v0.7.17-r1188
- Samtools 1.7
- Spades v3.13.0
- VarScan v2.4.4
- Pindel 0.2.5b9
- Delly 0.8.7
- snakemake 5.32.0

## Repo Contents

- meta-data 

  1. referece-genome，all the reference genomes used in the analysis
  2. assembly-genomes，genomes of P. aeruginosa isolates assembled by SPAdes.
  3. blastn results of recurrently mutated genes (blast to NCBI nt database)

- results:  result files for figure visualization and analysis 

- scripts: script files for analysis

## SNPs calling and annotation

1. Map the NGS reads to the reference genome (snakemake pipeline):

   `snakemake -s  scripts/analysis/NGS_map_PAO1-pipeline.py -p `

2. SNPs calling (snakemake pipeline):

   `snakemake -s scripts/analysis/iSNV_calling-pipeline.py -p ` 
   
   This script is used to call the SNPs and generate the summary file of the iSNVs and SNPs
   
3. Annotation of SNPs  by using SnpEff :

   - Convert the table obtained from "iSNV-calling" pipeline into vcf format(VCFv4.1) : 
     
     `python  scripts/analysis/iSNV_calling/iSNVTable_2_vcf_allsample.py  -i  $   -o $allsamplesVcfPath   -r   $ref_ID `
     
   - Annoting the iSNVs and SNPs :
     
      `java -jar snpEff.jar ann  $referenceID   $vcf_file  --protein  $ref_database_ID   >$outputFile_vcf`

## Repeat regions identification of the reference genome

- cut the ref genome to 150bp simulation sequences by step 1:

  `python  scripts/analysis/PA-refGnm_cutToSimmulateReads.py -r $refGnm -o $outputP -s $cutStep -l $readLength `

- align the simulated reads to the reference genome using blastn
1. build the database for blast
`makeblastdb  -in  $reference_genome_file   -input_type  fasta   -dbtype  nucl  -out   $databaseID`
2. Alignment using blastn
`blastn -db  $databaseID  -query  $simmulate_reads_fastaFile   -outfmt 6 -num_threads 4`
- Find the regions that may be repeat regions of the reference genome
  `python scripts/analysis/20210405-PA-simulateReads-mapToPArefAE004091-filt-RepeatRegion.py -i $blastOutfile  -o $output_file(PA.RepeatRegion-Cov30-identy60.Posi.txt) -r $referenceGenomeID `
  
- map the  simulation reads to the PAO1 reference genome 

  `snakemake -s scripts/analysis/simuReads_mapTo_PARefAE0049-pipeline.py -p`

## Identification of homologous regions with other bacteria in PAO1 genome 
- Cut the ref genome to 150bp simulate sequences:
  `python  scripts/analysis/PA-refGnm_cutToSimmulateReads.py -r $refGnm -o $outputP -s 1 -l 150 `
  
- Map the simulate sequences to the reference genomes of other bacteria

  `snakemake -s scripts/analysis/simuPA_mapTo_allBacteriaRefGnms-pipeline.py`

- Identify the homologous regions between the PAO1 reference genome of P.aeruginosa  and the reference genomes of other bacteria
`python scripts/analysis/PA-simulateReads-find-homoRegions.py -i $bamFile_simulateReadsMapToOtherbacteriaGenome -o $output_file -r $reference_ID`


## Phylogenetic analysis

- Phylogenetic analysis of samples based on the SNPs

`raxmlHPC-PTHREADS -s $aligned_seq_file -n raxmltree_result -m ASC_GTRGAMMA
--asc-corr=felsenstein -f a -x 12345 -N 1000 -p 123456 -T 24 -k`

- Phylogenetic analysis of core-genome sequences of samples and the public genomes

`fasttree -nt -gtr $aligned_sequence_file `



## The epidemic analysis of recurrent mutations in this study among public databases 

- Align the recurrently mutated gene sequences to the sequences of nt database by using NCBI blast program

- Filter the public sequences that have the same mutation with the recurrently mutated genes
scripts/analysis/find-the-Mutations-from -blastn-results.py

## Citation

If you use data, results or conclusion from this work, please cite:



## Acknowledgement



