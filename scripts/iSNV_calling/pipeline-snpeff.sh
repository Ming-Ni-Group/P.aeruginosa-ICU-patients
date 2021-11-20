#!/bin/bash

vcfP=$1
outSnpeffVcfP=$2
dbID=$3
snpeffSoftP=$4

#########################################################################3
echo 
echo "********************************************************"
echo "*                 mutation anno with snpeff            *"
echo "********************************************************"
echo

mkdir $outSnpeffVcfP

for vcf_F in $vcfP/*.vcf; do 
	echo "vcf_F:  "$vcf_F 
	fileID1=${vcf_F%.iSNVpy.vcf*}
	echo $fileID1
	fileID=${fileID1##*/}
	echo $fileID

	java -jar $snpeffSoftP/snpEff.jar  ann   -v  $dbID   $vcf_F > $outSnpeffVcfP/$fileID.snpeff.vcf

done


echo "================ done ===================="




