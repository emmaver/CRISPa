#!/bin/bash

#script for combining read files from separate lanes into one

#folder with raw sequencing data
path=$1
output_folder=$2


for folder in ` ls $path`;
do

#Concatenate multiple lanes into one file
cd "${path}${folder}"

for i in `ls *_R1_001.fastq.gz`;
do
	strain=`echo $i | rev | cut -d "/" -f1 | rev | cut -d "_" -f1`
	cat $i >> "${output_folder}/${strain}_R1.fastq.gz"
	echo "Done"
done


for i in `ls *_R2_001.fastq.gz`;
do
	strain=`echo $i | rev | cut -d "/" -f1 | rev | cut -d "_" -f1`
	cat $i >> "${output_folder}/${strain}_R2.fastq.gz"
	echo "Done"
done

done
