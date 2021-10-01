# Gene Per Locus Selection After Gene Scoring

Overall aim: to select the top gene  per locus post-GWAS 
- Genes have previously undergone machine learning, giving genes a ranking based on if they are likely to be causal for blood pressure
- Genes are grouped in loci (genomic regions  found to have statistical significance for blood pressure phenotypes) and so the aim of this code is to select the best gene within each group/locus based on outlined criteria below.

  Filtering to select gene(s) per locus if they are:
  1. RF_Score is > +1 standard deviation (SD). <br />
  If 1. is not met then:
  2. Genes > the average score are selected for further filtering by protein-protein interactions (PPI)
  3. Gene with highest direct PPI (direct_PPI_count) with known disease-genes is chosen OR if direct PPI is matching between genes they enter 4.
  4. Gene with highest secondary PPI (secondary_PPI_count) with known disease-genes is chosen OR if there are still matching PPI counts genes enter step 5.
  5. All genes matching PPIs are chosen

Example:

<img align="center" src="https://i.imgur.com/wiMDoaP.png" width="1000">
<br />
- Each locus is assigned a group number, e.g. the above genes are all in locus 616 so the selection is performed between those genes<br />

- No gene's machine learning score (RF_Score) is above the +1SD threshold<br />

- All genes > the average model score (MTHFR, NPPA, NPPB, and CLCN6) enter the PPI filtering steps, in which the genes also all match for their direct_PPI_count, but MTHFR has the largest secondary_PPI_count and so is the selected gene for that locus<br />

## How to run:
### Input files:
1. *Genelist_for_PPI.txt* - Gene list for identifying PPIs of genes of interest
2. *9606.protein.links.full.v11.0.txt* - PPI data needs to be downloaded from STRINGdb (https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens)
3. *BPLoci_ranked_Evangelou_RFR.txt* - Gene list with machine learning prioritisation scores per gene
4. *Total_genes_PPI_STRINGdb.txt* - Filtered PPIs from STRINGdb
	- PPI file is the output file of first script (STRINGdb PPI.r) and needed to run the main aims of the second script (Gene per locus selection.R)
	
### Files are input into two scripts<br />

**1. STRINGdb PPI.R:**  - getting PPIs for genes of interest <br />
- Files 1 (Genelist_for_PPI.txt) and 2 (9606.protein.links.full.v11.0.txt) enter STRINGdb PPI.R  - which aims to identify the STRING ids of input gene symbols and list the PPIs per gene, ultimately creating file 4 to run the second script

**2. Gene per locus selection.R**  - counting PPIs per gene and applying selection criteria <br />
- Files 3 (BPLoci_ranked_Evangelou_RFR.txt) and 4 (Total_genes_PPI_STRINGdb.txt) enter Gene per locus selection.R script which has 5 sections of code to run:
	1. Identify direct and secondary PPIs with known BP genes for all input genes
	2. Calculate statistics from model scores per locus (mean, SD, variance, upper SD etc.)
	3. Apply filtering rules to select genes per locus
	4. Plotting distribution of model scores for selected genes
	5. Count number of genes filtered at each step
	
## To Review:<br />
- Any and all feedback welcome<br />


```
> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252 
[2] LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] forcats_0.5.1     stringr_1.4.0     purrr_0.3.4       readr_1.4.0      
 [5] tibble_3.0.6      ggplot2_3.3.3     tidyverse_1.3.0   tidyr_1.1.2      
 [9] dplyr_1.0.4       data.table_1.13.6

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6        pillar_1.4.7      compiler_4.0.3    cellranger_1.1.0 
 [5] dbplyr_2.0.0      tools_4.0.3       jsonlite_1.7.2    lubridate_1.7.9.2
 [9] lifecycle_0.2.0   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10     
[13] reprex_1.0.0      cli_2.3.0         rstudioapi_0.13   DBI_1.1.1        
[17] haven_2.3.1       withr_2.4.1       xml2_1.3.2        httr_1.4.2       
[21] fs_1.5.0          generics_0.1.0    vctrs_0.3.6       hms_1.0.0        
[25] grid_4.0.3        tidyselect_1.1.0  glue_1.4.2        R6_2.5.0         
[29] readxl_1.3.1      modelr_0.1.8      magrittr_2.0.1    backports_1.2.1  
[33] scales_1.1.1      ellipsis_0.3.1    rvest_0.3.6       assertthat_0.2.1 
[37] colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.4      
[41] crayon_1.4.0   

```
