## Annotation contents

#### 1. acrfinder_adaptations

Contains the two Python scripts that were adapted from the original [AcrFinder tool](https://github.com/HaidYi/acrfinder), as well as the filtered database of trusted Aca genes.

#### 2. crisprcas_tool_comparison

Contains the Snakemake workflow for comparison of CASC, CRISPRCasFinder and CRISPRDetect for detection of CRISPR-Cas systems. A pdf file with snippets from their outputs and the Python scripts used to parse them are also provided.
     
#### 3. remaining files

The actual annotation pipeline (Snakefile) and its required auxiliary files, including an extended CRISPRCasFinder parser and two different AcrFinder parsers (a separate one for NCBI genomes).
