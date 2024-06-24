library(dplyr)
library(GEOquery)
library(tidyverse) 
library(DESeq2)
library(data.table)

##------GSE111016-----
## read in my counts table for the dataset
counts_GSE111016 <- fread('Count Data/GSE111016_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count111016 <- as.data.frame(counts_GSE111016, row.names = counts_GSE111016$GeneID)
gene_count111016 <- gene_count111016[, -1]

## series matrix file from geo
gse111016 <- getGEO('GSE111016', destdir = '.', getGPL = F)

## assesing needed metadata
meta_111016 <- pData(phenoData(gse111016[[1]]))

## sifting out my columns of interest
## manipulating variable to suit my downstream process using 'select' to choose my columns, 
#'rename' the columns to something easy to understand 
# mutate and gsub to rename column content to match count matrix 
## faced an error and had to use dplyr:: to specify which package will use rename
cond_111016<- meta_111016 %>%
  dplyr::select(19,11) %>%
  dplyr:: rename(sarcopenia = characteristics_ch1.1) %>%
  dplyr:: rename(samples = description) %>%
  mutate(samples = gsub('Sample ', 'Sample.', samples)) %>%
  mutate(sarcopenia = gsub('sarcopenia status: ', '', sarcopenia)) 

#matching my gene id to gene name 
## making rownames the sample id so it matches with my count colname
all(rownames(cond_111016) == colnames(gene_count111016))
all(colnames(gene_count111016) %in% rownames(cond_111016))

#setting factor levels
cond_111016$sarcopenia <- factor(cond_111016$sarcopenia)

## DDS object
dds_111016 <- DESeqDataSetFromMatrix(countData = gene_count111016, 
                              colData = cond_111016, design = ~ sarcopenia)

## prefiltering by removing rows with low gene count. 10 as the value of choice
dds_111016 <- dds_111016[names(dds_111016) %in% ahhesome_genes$entrezgene_id,] 

## variance stabilizing transformation
vsd_111016 <- vst(dds_111016, blind=FALSE)
vsd_111016_df <- as.data.frame(assay(vsd_111016))


##-----GSE111010-----
#GSE111010 DESEQ ANALYSIS 
count_GSE111010 <- fread('Count Data/GSE111010_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count_111010 <- as.data.frame(count_GSE111010, row.names = count_GSE111010$GeneID)
gene_count_111010 <- gene_count_111010[, -1]
## series matrix file from geo
gse111010 <- getGEO('GSE111010', destdir = '.', getGPL = F)

## assesing needed metadata
meta_111010 <- pData(phenoData(gse111010[[1]]))

## sifting out my columns of interest
## manipulating variable to suit my downstream process using 'select' to choose my columns, 
#'rename' the columns to something easy to understand 
# mutate and gsub to rename column content to match count matrix 
## faced an error and had to use dplyr:: to specify which package will use rename
cond_111010 <- meta_111010 %>%
  dplyr::select(11,12,13) %>%
  dplyr::rename(sarcopenia = characteristics_ch1.1) %>%
  dplyr::rename(low_muscle_mass = characteristics_ch1.2) %>%
  dplyr::rename(low_muscle_strength = characteristics_ch1.3) %>%
  mutate(low_muscle_mass = gsub('low muscle mass: ', '', low_muscle_mass)) %>%
  mutate(low_muscle_strength = gsub('low muscle strength and/or low physical performance: ', '', low_muscle_strength)) %>%
  mutate(sarcopenia = gsub('sarcopenia status: ', '', sarcopenia))

all(colnames(gene_count_111010) == rownames(cond_111010))
all(rownames(cond_111010) %in% colnames(gene_count_111010))

## set factor - multiple factors
cond_111010$sarcopenia <- factor(cond_111010$sarcopenia)
cond_111010$low_muscle_strength <- factor(cond_111010$low_muscle_strength) 

## create a deseq object to import count and sample conditions 
dds_111010 <- DESeqDataSetFromMatrix(countData = gene_count_111010, 
                                     colData = cond_111010, 
                                     design = ~ low_muscle_strength + sarcopenia)

## filtering with adhesome genes
dds_111010 <- dds_111010[names(dds_111010) %in% ahhesome_genes$entrezgene_id,] 

## variance stabilizing transformation
vsd_111010 <- vst(dds_111010, blind=FALSE) 
vsd_111010_df <- as.data.frame(assay(vsd_111010))


##----GSE111006----
#GSE111006 DESEQ ANALYSIS
count_GSE111006 <- fread('Count Data/GSE111006_raw_counts_GRCh38.p13_NCBI.tsv.gz') 
gene_count111006 <- as.data.frame(count_GSE111006, row.names = count_GSE111006$GeneID)
gene_count111006 <- gene_count111006[, -1]

## series matrix data
gse111006 <- getGEO('GSE111006', destdir = '.', getGPL = F) 

## assesing needed metadata
meta_111006 <- pData(phenoData(gse111006[[1]]))

## sifting out my columns of interest
## manipulating variable to suit my downstream process using 'select' to choose my columns, 
#'rename' the columns to something easy to understand 
# mutate and gsub to rename column content to match count matrix 
## faced an error and had to use dplyr:: to specify which package will use rename
cond_111006<- meta_111006 %>%
  dplyr::select(11,12,13) %>%
  dplyr:: rename(sarcopenia = characteristics_ch1.1) %>%
  dplyr:: rename(low_muscle_mass = characteristics_ch1.2) %>%
  dplyr:: rename(low_muscle_strength = characteristics_ch1.3) %>%
  mutate(low_muscle_mass = gsub('low muscle mass: ', '', low_muscle_mass)) %>%
  mutate(low_muscle_strength = gsub('low muscle strength and/or low physical performance: ', '', low_muscle_strength)) %>%
  mutate(sarcopenia = gsub('sarcopenia status: ', '', sarcopenia)) 


## matching rownames to colnames for deseqmatrix downstream -continued
all(rownames(cond_111006) == colnames(gene_count111006)) 
all(colnames(gene_count111006) %in% rownames(cond_111006))

## set factor - multiple factors
cond_111006$sarcopenia <- factor(cond_111006$sarcopenia)
cond_111006$low_muscle_strength <- factor(cond_111006$low_muscle_strength) 
cond_111006$low_muscle_mass <- factor(cond_111006$low_muscle_mass)

## create a deseq object to import count and sample conditions 
dds_111006 <- DESeqDataSetFromMatrix(countData = gene_count111006, 
                                     colData = cond_111006, 
                                     design = ~ low_muscle_strength +
                                       low_muscle_mass+ sarcopenia)

## prefiltering with adhesome genes
dds_111006 <- dds_111006[names(dds_111006) %in% ahhesome_genes$entrezgene_id,] 

## variance stabilizing transformation
vsd_111006 <- vst(dds_111006, blind=FALSE) 
vsd_111006_df <- as.data.frame(assay(vsd_111006))





##-----GSE174106----- 

count_GSE174106 <- fread('Count Data/GSE174106_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count174106 <- as.data.frame(count_GSE174106, row.names = count_GSE174106$GeneID)
gene_count174106 <- gene_count174106[, -1]

## series matrix import
gse174106 <- getGEO('GSE174106', destdir = '.', getGPL = F) 
## assesing needed metadata
meta_174106 <- pData(phenoData(gse174106[[1]])) 
## sifting out my columns of interest
## manipulating variable to suit my downstream process using 'select' to choose my columns,
cond_174106 <- meta_174106 %>%
 dplyr::select(13) %>%
  dplyr::rename(lean_mass = characteristics_ch1.3) %>% 
  mutate(lean_mass = gsub('appendicular lean mass (alm):', '', lean_mass, fixed = T)) %>%
  mutate(lean_mass = gsub ('stable low','stable_low', lean_mass)) %>%
  mutate(lean_mass = gsub ('stable high','stable_high', lean_mass))

all(rownames(cond_174106) %in% colnames(gene_count174106))
all(colnames(gene_count174106) == rownames(cond_174106))

## matching colnames of count data to rownames of series matrix 
## count column (samples) is less than coldata rownames 
## filter rownames by column to match colnmes to rownames
gene_count174106 <- subset(gene_count174106, select = -GSM5287346)
cond_174106 <- cond_174106[rownames(cond_174106) %in% colnames(gene_count174106), 
                                     ,  drop = F]

## setting my design factor
cond_174106$lean_mass <- factor(cond_174106$lean_mass)

## create a deseq object to import count and sample conditions 
dds_174106 <- DESeqDataSetFromMatrix(countData = gene_count174106, 
                                     colData = cond_174106, 
                                     design = ~ lean_mass) 


## filtering by removing rows with low gene count. 10 as the value of choice
dds_174106 <- dds_174106[names(dds_174106) %in% ahhesome_genes$entrezgene_id,] 

## variance stabilizing transformation
vsd_174106 <-  vst(dds_174106, blind=FALSE)
vsd_174106_df <- as.data.frame(assay(vsd_174106))



##----GSE159217------
## import counts and series matrix
count_GSE159217 <- fread('Count Data/GSE159217_raw_counts_GRCh38.p13_NCBI.tsv.gz') 
gene_count159217 <- as.data.frame(count_GSE159217, row.names = count_GSE159217$GeneID)
gene_count159217 <- subset(gene_count159217, select = -1) 

gse159217 <- getGEO('GSE159217', destdir = '.', getGPL = F)
meta_159217 <- pData(phenoData(gse159217[[1]]))

## manipulate metadata to create coldata 
cond_159217 <- meta_159217%>%
  dplyr::select(12)%>%
  dplyr::rename(age = characteristics_ch1.2) %>%
  mutate(age = gsub('age group: ','', age))

all(rownames(cond_159217) == colnames(gene_count159217))
all(colnames(gene_count159217) %in% rownames(cond_159217))

## set factor 
cond_159217$age <- factor(cond_159217$age)

##dds object 
dds_159217 <- DESeqDataSetFromMatrix(countData = gene_count159217, 
                       colData = cond_159217, design = ~ age)
## filter with adhesome genes
dds_159217 <- dds_159217[names(dds_159217) %in% ahhesome_genes$entrezgene_id,] 

## variance 
vsd_159217 <- vst(dds_159217, blind = F)
vsd_159217_df <- as.data.frame(assay(vsd_159217))

 



##-----GSE163434----

## import count and series matrix
count_GSE163434 <- fread('Count Data/GSE163434_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count163434 <- as.data.frame(count_GSE163434, row.names = count_GSE163434$GeneID)
gene_count163434 <- gene_count163434[, -1]

gse163434 <- getGEO('GSE163434', destdir = '.', getGPL = F)
meta_163434 <- pData(phenoData(gse163434[[1]]))

cond_163434 <- meta_163434 %>%
  dplyr::select(12)%>%
  dplyr::rename(time_point = characteristics_ch1.2) %>%
  mutate(time_point = gsub('time point: ','',time_point))

## filtering out the rows post 14 weeks
gene_count163434 <- gene_count163434[, -c(2,4)]
## filter rownames by colnames
cond_163434 <- cond_163434[rownames(cond_163434)  %in% colnames(gene_count163434), , drop = F]

dds_163434 <- DESeqDataSetFromMatrix(countData = gene_count163434,
                                     colData = cond_163434, design = ~ 1)
## filter by adhesome genes 
dds_163434 <- dds_163434[names(dds_163434) %in% ahhesome_genes$entrezgene_id,] 

## varaince stabilization transformation
vsd_163434 <- vst(dds_163434, blind = F) 
vsd_163434_df <- as.data.frame(assay(vsd_163434))

