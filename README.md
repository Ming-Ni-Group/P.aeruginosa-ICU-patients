# title



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

- meta-data: sample data and preprocessed data used for analysis.

  1. PA ref genome
  2. genomes of the other bacteria 
  3. samples' information file

- results:  The results from the 

- scripts: for analysis and figures plot

## iSNVs and SNPs calling and annotation

1. Map the WGS reads to the reference genome :

   using snakemake 

   `snakemake -s  NGS_map_PAO1-pipeline.py`

2. iSNVs and SNPs calling :

   `snakemake -s iSNV_calling-pipeline.py `
   
   (This can also generate the stat of the iSNVs and SNPs, samples' genomes replaced the SNPs cites)
   
3. Annotation of iSNVs and SNPs  by using SnpEff :

   - Convert iSNV table into vcf format(VCFv4.1) : 
     
     `python  iSNVTable_2_vcf_allsample.py  -i  $   -o $allsamplesVcfPath   -r   MN908947.3`
     
   - Annoting the iSNVs and SNPs :
     
      `java -jar snpEff.jar ann  $referenceID   $vcf_file    >$outputFile_vcf`

## Repeat regions in the reference genome

- cut the ref genome to 150bp simulation sequences:

  `python  PA-refGnm_cutToSimmulateReads.py -r $refGnm -o $outputP -s $cutStep -l $readLength `

- map the  simulation reads to the other bacteria genomes

  `snakemake -s simuReads_mapTo_PARefAE0049-pipeline.py`

- align the simulated reads to the reference genome using blastn
1. build the database for blast
`makeblastdb  -in  $reference_genome_file   -input_type  fasta   -dbtype  nucl  -out   $databaseID`
2. Alignment using blastn
`blastn -db  $databaseID  -query  $simmulate_reads_fastaFile   -outfmt 6 -num_threads 4`
- Find the regions that may be repeat regions of the reference genome
  `python 20210405-PA-simulateReads-mapToPArefAE004091-filt-RepeatRegion.py -i $blastOutfile  -o $output_file(PA.RepeatRegion-Cov30-identy60.Posi.txt) -r $referenceGenomeID `

## Identity the homo regions with other bacteria
- Cut the ref genome to 150bp simulate sequences:
  `python  PA-refGnm_cutToSimmulateReads.py -r $refGnm -o $outputP -s 1 -l 150 `
  
- Map the simulate sequences to the reference genomes of other bacteria

  `snakemake -s simuPA_mapTo_allBacteriaRefGnms-pipeline.py`

- Find the homo regions with other bacteria
`python PA-simulateReads-find-homoRegions.py -i $bamFile_simulateReadsMapToOtherbacteriaGenome -o $output_file -r $reference_ID`


## Phylogenetic analysis

- Prepare the  aligned sequences based on SNPs cites:

` python generate-sequences-of SNPs-cites.py -i $iSNV_SNP_table_Path  -f $Samps_info_file    -s  $path_of_Annotionedfiles_by_snpeff  -R $repeatRegionsCites -E $HomoRegionsCites  -o $output_path    -M  0.9  -n $ntfreqFiles_path`

- Phylogenetic analysis

`raxmlHPC-PTHREADS -s $align_SNPs_file -n raxmltree_result -m GTRGAMMAI -f a -x 12345 -N 1000 -p 123456 -T 24 -k`

## Genome assembly

genome assembly was conducted by using Spades:

`spades.py -o  $output_path  -t 30  -k 21,33,45,55,67,77,89,99,103,111,115,127  --careful  -1  $reads_R1   -2   $reads_R2`



## Identity of Indels 

- Indels calling with Pindel:

  `snakemake -s `

- Indels calling with Delly

- Filt the Indels (shared by both results of the two softwares)



## Citation

If you use data, results or conclusion from this work, please cite:



## Acknowledgement



