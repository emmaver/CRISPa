## Assembly contents

#### 1. ncbi_genomes

Contains the list of RefSeq FTP addresses of the genomes that passed the GC % and genome size filters, as well as a bash script which was used to fetch the assembly and annotation files from these genomes.

#### 2. short_read_assembler_comparison

Contains the Snakemake workflow for comparison of SPAdes (default), SPAdes (--careful) and Shovill for assembly of short Illumina reads.
     
#### 3. remaining files

`concatenate_sequencing_lanes.sh` is a script used for consolidation of the sequencing reads over multiple lanes into one *fastq* file. The other files comprise the actual assembly pipeline (Snakefile) and its required auxiliary files.
