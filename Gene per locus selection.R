#################### Gene Per Locus Selection After Gene Priotiziation ####################
#*Restart R after running STRINGdb PPI.R if running this code straight after

#1. Identify direct and secondary PPIs with known BP genes for all genes
#2. Calculate statistics from model scores per locus (mean, SD, variance, SD+1 etc.)
#3. Apply filtering rules to select genes per locus
#4. Plotting distribution of model scores for selected genes
#5. Count number of genes filtered at each step


#setwd("~/PhD Year 2/")
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

#####################################################################################
#1. Identify direct and secondary PPIs with known BP genes for all input scored genes

#Read in PPI counts from interactions in stringdb
interactors <- fread("Total_genes_PPI_STRINGdb.txt") #file is made from downloadable STRINGdb PPI data which has been processed to list interactors per gene

#Read in file of genes grouped by loci
bploci <- fread("BPLoci_ranked_Evangelou_RFR.txt")
colnames(bploci)[4] <- 'RF_Score'

#Merge files to give each BP loci gene its matching interactors listed in an additonal column
bpmerged  <- merge(bploci, interactors, by='Gene', all.x=TRUE)

#Identify known BP genes - to later find their interactors (to be used in further gene filtering below)
bpknown <- filter(bpmerged, Training_Score == 1)

#Identify genes with a direct interaction with known BP genes (creating 'direct_PPI_count' column)   
bp_PPI1 <- bpmerged %>%
  mutate(direct_PPI_count = str_count(interactors, str_c(bpknown$Gene, collapse = '|')))

#Genes with no STRINGdb identifiers have NA for PPI counts and need to be set to 0
bp_PPI1$direct_PPI_count[is.na(bp_PPI1$direct_PPI_count)] <- 0

#Identify genes that have interactors interacting with interactors of BP genes (creating 'secondary_PPI_count' column)

#Interactors need to be separated to be individually counted (initially all interactors are comma separated in one cell per gene)
sep_rows <- . %>% separate_rows(interactors, sep = ", ")

#Code for bp_PPI2 will produce an empty object if ran straight after STRINGdb PPI code - need to restart R
bp_PPI2 <- bp_PPI1 %>% 
  sep_rows() %>% 
  mutate(
    found = !is.na(match(interactors, sep_rows(bpknown)$interactors))) %>% 
  group_by(Gene) %>% 
  summarise(
    interactors = toString(interactors), 
    secondary_PPI_count = sum(found)
  )

#Combining direct and secondary PPI counts into 1 dataset and dropping interactors column
bp_PPI_total <- merge(bp_PPI1, bp_PPI2,by='Gene', all.x=TRUE)
bp_PPI_total <- select(bp_PPI_total, Gene, direct_PPI_count, secondary_PPI_count)

#Getting genes at their loci with their PPI counts
merged <- merge(bploci, bp_PPI_total, by='Gene', all.x=TRUE)
merged <- unique(merged)

#Re-ordering columns to more easily view genes and their PPI counts
df <- select(merged, loci, Gene, RF_Score, Training_Score, direct_PPI_count, secondary_PPI_count, LocusName, rsID_LEAD_SNP, CP_LEADSNP, CPs,
             Source, Type,  min_pvalue, GPrior_Score, ToppGene_Score, Mantis_probability, minTRAIT, GTEx_signif_sexbias)

############################################################################################################
#2. Calculate statistics from model scores per locus (to get +1SD - Upper_SD_Threshold' column for filtering)

#Creating statistical columns from model predictions
df <- df %>% group_by(loci) %>% mutate(AvgScore_Per_Locus=mean(RF_Score,na.rm = T),
                                       diff = max(RF_Score) - RF_Score,
                                       Avg_TopGene_Diff_Per_Locus = sum(diff) / (n() - 1),
                                       Variance_Per_Locus = var(RF_Score,na.rm = T),
                                       IQR_Per_Locus = IQR(RF_Score,na.rm = T),
                                       SD_Per_Locus = sd(RF_Score,na.rm = T),
                                       LocusGenes = toString(Gene))

df$Upper_SD_Threshold <- df$AvgScore_Per_Locus + df$SD_Per_Locus

###################################################
#3. Apply filtering rules to select genes per locus

  #Filtering to select gene(s) per locus if they are:
  #1. >Upper_SD Threshold (+1 SD) OR
  #2. If multiple genes between average score - Upper SD then filter by PPI
  #3. Gene with highest direct PPI chosen OR if matching direct_PPI_count:
  #4. Gene with highest secondary PPI chosen OR if matching secondary_PPI_count:
  #5. all genes left chosen


#Identifying genes that meet the first condition of being > +1SD (Upper_SD_Threshold):
UpperSD_Genes <- df %>%
  group_by(loci) %>%
  #filtering within loci with multiple genes (loci with already singular genes not entering filter):
  filter(if(n() > 1) {(RF_Score > Upper_SD_Threshold) } else TRUE) %>%
  filter(if(n() > 1) {RF_Score > Upper_SD_Threshold & direct_PPI_count == max(direct_PPI_count)} else TRUE) %>%
  filter(if(n() > 1) {RF_Score > Upper_SD_Threshold & secondary_PPI_count == max(secondary_PPI_count)} else TRUE)

#Identifying loci where only 1 gene is > average model score (and therefore no further PPI filtering is needed):
Avg_Genes <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  filter(RF_Score > AvgScore_Per_Locus) %>%
  filter(n() == 1)

#Selecting gene(s) per locus depending direct and secondary PPI filters:
PPI_Genes <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  subset(!(loci %in% Avg_Genes$loci)) %>%
  filter(RF_Score >= AvgScore_Per_Locus) %>%
  filter(direct_PPI_count == max(direct_PPI_count)) %>%
  filter(secondary_PPI_count == max(secondary_PPI_count)) 

#Combining selected genes per loci:
new_df <- rbind(UpperSD_Genes, Avg_Genes, PPI_Genes)
new_df <- new_df %>% arrange(loci)

final <- select(new_df, loci, Gene, RF_Score, Training_Score, LocusGenes, AvgScore_Per_Locus,
                SD_Per_Locus,Upper_SD_Threshold, direct_PPI_count, secondary_PPI_count,
                Avg_TopGene_Diff_Per_Locus,
                Variance_Per_Locus, IQR_Per_Locus, LocusName, rsID_LEAD_SNP,
                Source, Type, CP_LEADSNP, CPs, minTRAIT, min_pvalue, GPrior_Score, ToppGene_Score, Mantis_probability, GTEx_signif_sexbias)

colnames(final)[1] <- "Loci"
colnames(final)[3] <- "RF Score"
colnames(final)[5] <- "LocusGene(s)"
final <- unique(final)

write.csv(final, "./Gene_per_locus_PPI_SD_Evangelou_RF.csv", row.names = FALSE)

#########################################################
#4. Plotting selected genes per locus score distribution:

ggplot(final, aes(x=`RF Score`))+
  geom_histogram(color="darkblue", fill="lightblue") +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggtitle("RF Prioritisation Score Distribution (Selected Genes Per Locus)")  +  
  ylab("Gene Count") +
  theme(text = element_text(size=16), plot.title = element_text(size = 16, hjust = 0.5))

ggsave("BPLoci_Distribution_selected_Evangelou_RF.tiff",width = 10, height = 4, dpi=300,
       limitsize = FALSE)

###############################################
#5. Count number of genes filtered at each step

#Counting number of loci with >1 gene selected after all filtering steps:
multigeneloci <- new_df[duplicated(new_df$loci) | duplicated(new_df$loci, fromLast = TRUE), ] 
length(unique(multigeneloci$loci)) 
#write.csv(multigeneloci , "./duploci.csv", row.names = FALSE) 

#Number of genes that were selected by 1st filter (> +1 SD model score) are in UpperSD_Genes object
nrow(UpperSD_Genes)
#Number of genes that were selected by being the only gene > the average score are in Avg_Genes object
nrow(Avg_Genes)

#Counting number of loci entering PPI filtering:
length(unique(PPI_Genes$loci)) 

#Count number of loci (seeing how many still have >1 gene past direct PPI filter and are entering final secondary PPI filter)
count1stPPI <- df %>% 
  group_by(loci) %>% 
  subset(!(loci %in% UpperSD_Genes$loci)) %>%
  subset(!(loci %in% Avg_Genes$loci)) %>%
  filter(RF_Score >= AvgScore_Per_Locus) %>%
  filter(direct_PPI_count == max(direct_PPI_count))

count1stPPI1 <- select(count1stPPI, Gene, loci)
#Counting number of loci entering 2nd PPI filter:
PPI1_dup <- unique(count1stPPI1[duplicated(count1stPPI1$loci), "loci"])  
#Counting number of loci entering selected by 1st PPI filter:
count1stPPI2 <-count1stPPI1[!count1stPPI1$loci %in% PPI1_dup$loci,] 

