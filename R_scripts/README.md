This folder contains 4 scripts which follow each other logically (the dataset is gradually built up).

#### 1. asm_data_consolidation.R

First script in which the NCBI genomes are filtered, genome size distribution is constructed and assembly results are assessed.

#### 2. CRISPR-Cas.R

Second script in which CRISPR-Cas annotations are explored and size distributions reconstructed. Phylogenetic group and host organism annotations are consolidated as well.

#### 3. anti-CRISPR.R

Third script in which anti-CRISPR results are explored and combined with CRISPR-Cas to obtain *corrected* annotations.
Size distributions etc. are again considered.

#### 4. infectivity.R

Final script in which the CRISPR-Cas and anti-CRISPR annotations are assessed in light of phage infectivity data.

---

### List of packages used over all scripts

* caret
* dplyr
* ggboxplot
* ggfortify
* ggplot2
* ggpubr
* gridExtra
* plyr
* readxl
* tree
