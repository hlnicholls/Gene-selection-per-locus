####Protein-Protein Interactions for GWAS associated BP genes####
#1. Match genes to their STRINGdb identifiers
#2. Match gene and string identifiers to interactors per gene 
#3. Collapsing interactors to be listed per gene in each row

#setwd("~/PhD Year 2/")
library(data.table)
library(STRINGdb)
library(dplyr)
library(plyr)
library(tidyverse)
library(sqldf)

#############################################
#1. Match genes to their STRINGdb identifiers

Genelist <- fread("Genelist_forPPI.txt")

#Using STRINGdb's R package set to newest string version (11) and human species by NCBI taxonomy number (9606)
#Set score_threshold to 0 to view all possible interactions (score is string's combined score from all their sources for PPIs)
string_db <- STRINGdb$new(version="11", species=9606,
                          score_threshold=0, input_directory="")

#Identify gene's STRINGdb identifiers:
#string_db$map function adds an additional column with STRING identifiers to the dataframe that is passed as first parameter (gene list)
genes_mapped <- string_db$map(Genelist, "Gene", removeUnmappedRows = TRUE)

#genes_mapped
#Gene   STRING_id
#A1CF	  9606.ENSP00000378868
#A4GALT	9606.ENSP00000384794
#AACS	  9606.ENSP00000324842


#############################################################
#2. Match gene and string identifiers to interactors per gene 

#Total data downloaded from STRING including all score types per interaction loaded
#Data in proteindf and was downloaded from: https://string-db.org/cgi/download?sessionId=bLbfFhAKfphK

proteindf <- fread("9606.protein.links.full.v11.0.txt")
dim(proteindf)
#[1] 11759454       16
names(proteindf)
#[1] "protein1"                 "protein2"                 "neighborhood"            
#[4] "neighborhood_transferred" "fusion"                   "cooccurence"             
#[7] "homology"                 "coexpression"             "coexpression_transferred"
#[10] "experiments"              "experiments_transferred"  "database"                
#[13] "database_transferred"     "textmining"               "textmining_transferred"  
#[16] "combined_score" 

colnames(proteindf)[1] <- "STRING_id1"
colnames(proteindf)[2] <- "STRING_id2"

#Interaction scores only from coexpression, experiments and database columns used:
protdf <- select(proteindf, STRING_id1, STRING_id2, coexpression, experiments, database)

#Duplicating STRING_id in genes_mapped object to have renamed matching column STRING_id1 
  #- Aim is to match with named STRING_id1 interactors in STRINGdb dataset, finding the gene names for those interactors
genes_mapped$STRING_id1 <- genes_mapped$STRING_id 
genes_mapped$STRING_id2 <- genes_mapped$STRING_id 

#Identifying gene names for first STRING protein ID column:
mapped1 <- select(genes_mapped, STRING_id1, Gene)
colnames(mapped1)[2]<- "Gene1"
string1df <- merge(protdf, mapped1, by='STRING_id1', all.y=TRUE)
colnames(string1df)[6]<- "Gene1"
#head(string1df)
#STRING_id1           STRING_id2          coexpression experiments database  Gene1
#9606.ENSP00000002596 9606.ENSP00000306425            0           0        0 HS3ST1
#9606.ENSP00000002596 9606.ENSP00000261523            0           0        0 HS3ST1
#9606.ENSP00000002596 9606.ENSP00000347088            0           0        0 HS3ST1
#9606.ENSP00000002596 9606.ENSP00000389176            0           0        0 HS3ST1
#9606.ENSP00000002596 9606.ENSP00000216124            0           0        0 HS3ST1
#9606.ENSP00000002596 9606.ENSP00000199447            0           0        0 HS3ST1


#Identifying gene names for second STRING protein ID column:
mapped2 <- select(genes_mapped, STRING_id2, Gene)
string2df <- merge(protdf, mapped2, by='STRING_id2', all.y=TRUE)
colnames(string2df)[6]<- "Gene2"

mergedIDs <- inner_join(string1df, string2df)
mergedIDs <- select(mergedIDs, STRING_id1, Gene1, STRING_id2, Gene2, coexpression, experiments, database)
#mergedIDs has matched gene symbols for string IDs in both interactor columns of the orignal protein data downloaded from string:
#dim(mergedIDs)
#2331566       7
#head(mergedIDs)
#STRING_id1           Gene1           STRING_id2   Gene2 coexpression experiments database
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000306425  CHCHD7            0           0        0
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000261523    RORA            0           0        0
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000347088   LARGE            0           0        0
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000389176 C1GALT1            0           0        0
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000199447    NME8            0           0        0
#9606.ENSP00000002596 HS3ST1 9606.ENSP00000264893  SEPT11            0           0        0


#Filtering interactions by STRING's lowest confidence threshold and 3 measurements:
df2filt <- filter(mergedIDs , experiments!= 0 | coexpression!=0 | database!=0 )
df2filt$Max <- apply(df2filt[,5:7], 1, FUN=max)
d2filt <- filter(df2filt, Max >= 150)

############################################################
#3. Collapsing interactors to be listed per gene in each row
# Duplicated genes in above output examples needs to be condensed to 1 row with all counterpart interactors listed for it 

#Compress gene's interactors for first Gene1 column:
genes1 <- df2filt[, c(2,4)]
genes1<-genes1[!(genes1$Gene1=="NA"),]
genes1<- data.table(genes1)
genes1 <- genes1[,lapply(.SD, function(col) paste(col, collapse=", ")), 
                 by=.(Gene1)]
genes1$Gene2 <- gsub("NA,", "", genes1$Gene2)
#First row of genes1:
#Gene1  #Gene2
#HS3ST1	CACNG2, EHF, GPC6, SDC2, ...

#Compress gene's interactors for second Gene2 column:
genes2 <- df2filt[, c(2,4)]
genes2<-genes2[!(genes2$Gene2=="NA"),]
genes2<- data.table(genes2)
genes2 <- genes2[,lapply(.SD, function(col) paste(col, collapse=", ")), 
                 by=.(Gene2)]
genes2$Gene1 <- gsub("NA,", "", genes2$Gene1)
#First row of genes2:
#Gene2  Gene1
#CACNG2	HS3ST1, CACNG3, LRRC7, ...

#Combine into one dataset:
colnames(genes1)[1] <- "Gene"
colnames(genes2)[1] <- "Gene"

interactdf <- merge(genes1, genes2, by='Gene', all=T)
#Combine interactors from both gene1 and gene2 and remove duplicates:
interactdf  <- transform(interactdf, interactors = paste(Gene1, Gene2, sep = ","))
interactdf$interactors <- sapply(interactdf$interactors , function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

interactors <- select(interactdf, Gene, interactors)

write.table(interactors,"Total_genes_PPI_STRINGdb.txt",sep="\t",row.names=F, quote=F)

