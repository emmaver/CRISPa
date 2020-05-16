#!/bin/bash
WD=/scratch/leuven/325/vsc32567/MT/crispr_annotation/NCBI/genomes

REFSEQ=$WD/../pass_genomes_refseq.txt

B=_assembly_stats.txt
C=_genomic.fna.gz
D=_genomic.gff.gz
E=_genomic.gbff.gz
J=_protein.faa.gz
G=_feature_count.txt.gz
H=_feature_table.txt.gz

while read FTP
do
    F=$(tr -dc '[[:print:]]' <<< "$FTP")
    ISOLATE=$(basename $F)
    I=$(tr -dc '[[:print:]]' <<< "$ISOLATE")
    DIR=$WD/$I
    mkdir -p $DIR
    cd $DIR
    ADDR=$F/$I
    A=$(tr -dc '[[:print:]]' <<< "$ADDR")

    wget --quiet --timestamping "$A$B"
    wget --quiet --timestamping "$A$C"
    wget --quiet --timestamping "$A$D"
    wget --quiet --timestamping "$A$E"
    wget --quiet --timestamping "$A$J"
    wget --quiet --timestamping "$A$G"
    wget --quiet --timestamping "$A$H"

    echo "Saved isolate $I in folder $DIR"

done < $REFSEQ
