Below is a list of remaining commands that were used in the project (outside of the Snakemake workflows)

## Kraken analysis
To test some suspicious isolates for contamination, the taxonomic classification tool Kraken 2 was used.

First, a custom database of all bacterial, viral and plasmid sequences was built: 

```
kraken2-2.0.8-beta/kraken2-build --download-taxonomy --db bacvirplasm

kraken2-2.0.8-beta/kraken2-build --download-library bacteria --db bacvirplasm 
kraken2-2.0.8-beta/kraken2-build --download-library plasmid --db bacvirplasm 
kraken2-2.0.8-beta/kraken2-build --download-library viral --db bacvirplasm 

kraken2-2.0.8-beta/kraken2-build --build --db bacvirplasm 
```

Next, this database was used to classify the supplied assemblies (based on k-mer sequences):

```
kraken2 --db bacvirplasm assembly.fasta --output file --report file
```


## Roary pan-genome analysis

Two analyses were performed with the following command (analysis A for a subset of 603 genomes, analysis B for all 5382 isolates):

```
roary -p 36 -g 60000 -f ./roary_5382 -e -n -v ../crispr_annotation/annotation/all_gff/*gff
```

## FastTree phylogenetic tree construction

The core gene alignment produced by Roary for analysis A was supplied to FastTree to infer the population structure of these 603 genomes.
FastTreeMP is the multithreaded version of FastTree.

```
./FastTreeMP -nt -gtr -log ft_fast.log -nosupport -fastest < core_gene_alignment_603.aln > 603_fast_tree.nw
```
