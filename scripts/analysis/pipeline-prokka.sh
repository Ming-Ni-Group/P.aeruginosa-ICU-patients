#!/bin/bash

gffP=$1
outP=$2
genbankdb=$3  ## tonglv --- /home/amax/anaconda3/envs/biotools_lhj/db/Pseudomonas/AE004091.2.gb


#########################################################################3
echo 
echo "********************************************************"
echo "*                 prokka genome anno                   *"
echo "********************************************************"
echo


for gff_F in $gffP/*.fasta; do 
	echo "gff_F:  "$gff_F 
	ID1=${gff_F##*/}
	ID=${ID1%.fasta*}

	prokka  --kingdom Bacteria  --force  --outdir  $outP/${ID}    --prefix  $ID  \
	--proteins   $genbankdb      $gff_F
done


echo "================ done ===================="




