getwd()
library(dplyr)
library(GEOquery)
library(tidyverse) 
library(DESeq2)
library(data.table)
library(WGCNA)
library(CorLevelPlot)

## reading in the specific muscle adhesome gene metadata
ahhesome_genes <- read.csv('Meta_Adhesome.csv')

##--------GSE226151--------------------------------
##read my dataset into R
count_GSE226151 <- fread('Count Data/GSE226151_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count226151 <- as.data.frame(count_GSE226151, row.names = count_GSE226151$GeneID)
gene_count226151 <- gene_count226151[,-1]

## read the series matrix using geoquery
gse226151 <- getGEO('GSE226151', destdir = '.', getGPL = F)

## got the metadata containing file info
meta_226151 <- pData(phenoData(gse226151[[1]]))

## sifting out my columns of interest
## manipulating variable to suit my downstream process using 'select' to choose my columns, 
#'rename' the columns
# mutate and gsub to rename column content to match count matrix
cond_226151 <- meta_226151 %>%
  dplyr::select(1,13,14) %>% 
  dplyr::rename(disease = title) %>%
  dplyr::rename(age = characteristics_ch1.4) %>%
  dplyr::rename(status = characteristics_ch1.3) %>%
  mutate(disease = gsub('Skeletal muscle,', '', disease)) %>%
  mutate(age = gsub('age:', '', age)) %>%
  mutate(disease = gsub(' 0', '_', disease)) %>%
  mutate(disease = gsub('HA ', 'HA_', disease)) %>% 
  mutate(disease = gsub('PS ', 'PS_', disease)) %>%
  mutate(disease = gsub('S ', 'S_', disease)) %>%
  mutate(status = gsub('disease state:', '', status)) %>%
  mutate(status = gsub('healthy aged', 'healthy_aged', status)) %>%
  mutate(status = gsub('pre-sarcopenia', 'pre_sarcopenia', status))


## filter rownames of cond df by colnames of gene count
cond_226151 <- cond_226151[rownames(cond_226151) %in% colnames(gene_count226151), , drop = F]

## confirming if rowname matches colname. first code aligns all row and col names
all(colnames(gene_count226151) %in% rownames(cond_226151)) 
all(colnames(gene_count226151) == rownames(cond_226151))

## set factor
cond_226151$status <- factor(cond_226151$status)

##creating a Deseqdataset
dds_226151 <- DESeqDataSetFromMatrix(countData = gene_count226151, 
                              colData = cond_226151, design = ~ 1 + status)

## filtering dds by adheosme genes
dds_226151 <- dds_226151[names(dds_226151) %in% ahhesome_genes$entrezgene_id,]

## variance stabilizing transformation
vsd_226151 <- vst(dds_226151, blind = FALSE) 
vsd_226151_df <- as.data.frame(assay(vsd_226151))
 

 
 
 
 
 


##----------GSE242202------------ 
## DDS GSE242202
count_GSE242202 <- fread('Count Data/GSE242202_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count242202 <- as.data.frame(count_GSE242202, row.names = count_GSE242202$GeneID)
gene_count242202 <- gene_count242202[,-1]

## reading in my series matrix 
gse242202 <- getGEO('GSE242202', destdir = '.', getGPL = F) 

## got the metadata containing file info
meta_242202 <- pData(phenoData(gse242202[[1]])) 

## sifting out my columns of interest
## manipulating variable to suit my downstream process conditions data 
cond_242202 <- meta_242202 %>%
  dplyr::select(1,12,13) %>% 
  dplyr::rename(age = characteristics_ch1.2) %>%
  dplyr::rename(disease = characteristics_ch1.3) %>%
  mutate(age = gsub('age: ', '', age)) %>%
  mutate(disease = gsub('healthy/patient: ', '', disease))

## ensuring rownames and colnames align - some columns missing in counts, so I filtered to remove them
cond_242202 <- cond_242202[rownames(cond_242202)  %in% colnames(gene_count242202), , drop = F]

all(colnames(gene_count242202) == rownames(cond_242202))
all(rownames(cond_242202)  %in% colnames(gene_count242202))
all(rownames(cond_242202)  == colnames(gene_count242202))

## create dds object 
dds_242202 <- DESeqDataSetFromMatrix(countData = gene_count242202, 
                       colData = cond_242202,
                       design = ~ disease) 

## filtering by adhesome genes
dds_242202 <- dds_242202[names(dds_242202) %in% ahhesome_genes$entrezgene_id,]

## variance stabilizing transformation
vsd_242202 <- vst(dds_242202, blind = F)  
vsd_242202_df <- as.data.frame(assay(vsd_242202))





##---------GSE235781-----
## DDS GSE235781
count_GSE235781 <- fread('Count Data/GSE235781_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count235781 <- as.data.frame(count_GSE235781, row.names = count_GSE235781$GeneID)
gene_count235781 <- gene_count235781[,-1]

## series matrix import
gse235781 <-getGEO('GSE235781', destdir = '.', getGPL = F)
meta_235781 <- pData(phenoData(gse235781[[1]])) 

## sifting out my columns of interest form series matrix 
cond_235781 <- meta_235781 %>%
  dplyr::select(1,11,12,13,16) %>%
  dplyr::rename(fiber_type = characteristics_ch1.1)%>%
  dplyr:: rename(age = characteristics_ch1.2) %>%
  dplyr::rename(training = characteristics_ch1.3) %>%
  dplyr::rename(time_point = characteristics_ch1.6) %>%
  dplyr::mutate(fiber_type = gsub('fiber-type: ', '', fiber_type)) %>%
  mutate(age = gsub('age.group: ', '', age)) %>%
  mutate(training = gsub('training.state: ', '', training)) %>%
  mutate(time_point = gsub('time.point: ', '', time_point))
 
## colnames of counts df was shorter than rownames of series matrix data 
## filtered rownames with colnames to ensure they match 
cond_235781 <- cond_235781[rownames(cond_235781)  %in% colnames(gene_count235781), , drop = F]
##deleting unneeded data from series matrix - everything that has 'POST workout' in rows 
cond_235781 <- cond_235781[-which(cond_235781$time_point == 'POST'),] 

colnames(gene_count235781)  %in% rownames(cond_235781) 
colnames(gene_count235781)  == rownames(cond_235781)
rownames(cond_235781) %in% colnames(gene_count235781)

## identify columns that are present as rownames in cond df - and filter count df by it
gene_count235781<- gene_count235781[, colnames(gene_count235781)  %in% rownames(cond_235781)]

## manipulating data be understandable
cond_235781 <- cond_235781 %>%
  mutate(age = gsub('LLE', 'old', age)) %>%
  mutate(age = gsub('YE', 'young', age))

## factor - age
cond_235781$age <- factor(cond_235781$age)

## deseq object - used round cus I faced an error - some values in assay are not integers
dds_235781 <- DESeqDataSetFromMatrix(countData = round(gene_count235781),
                       colData = cond_235781,
                       design = ~ age) 
## filtering deseq object with adhesome genes 
dds_235781 <- dds_235781[names(dds_235781) %in% ahhesome_genes$entrezgene_id,] 

##variance stabilization transformation
vsd_235781 <- vst(dds_235781, blind = F)
vsd_235781_df <- as.data.frame(assay(vsd_235781))








##----------GSE144304----------
## DDS FOR GSE144304
count_GSE144304 <- fread('Count Data/GSE144304_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count144304 <- as.data.frame(count_GSE144304, row.names = count_GSE144304$GeneID)
gene_count144304 <- gene_count144304[, -1]

## series matrix import from geo
gse144304 <- getGEO('GSE144304', destdir = '.', getGPL = F) 
meta_144304 <- pData(phenoData(gse144304[[1]]))

## series data manipulation to make a conditions table 
cond_144304 <- meta_144304 %>%
  dplyr::select(1,11,13) %>%
  dplyr::rename(age = characteristics_ch1.1) %>%
  dplyr::rename(muscle_status = characteristics_ch1.3) %>%
  mutate(age = gsub('age: ', '', age)) %>%
  mutate(muscle_status = gsub('group: ', '', muscle_status)) %>%
  mutate(age = gsub('\\..*', '', age))

colnames(gene_count144304) %in% rownames(cond_144304)
rownames(cond_144304) %in% colnames(gene_count144304)
all(colnames(gene_count144304) == rownames(cond_144304))

## setting factor and levels for 
cond_144304$muscle_status <- factor(cond_144304$muscle_status, levels = c('Young', 
                                                                        'Frail', 'Fit')) 


## make dds object - design set for 3 levels
dds_144304 <- DESeqDataSetFromMatrix(countData = gene_count144304, 
                       colData =cond_144304,
                       design = ~ muscle_status)

#filtering ddds 
dds_144304 <- dds_144304[names(dds_144304) %in% ahhesome_genes$entrezgene_id,]

## variance stabilizing transformation
vsd_144304 <-  vst(dds_144304, blind=FALSE) 
vsd_144304_df <- as.data.frame(assay(vsd_144304))






##------------GSE167186- no count txt------------
## DDS for GSE167186
count_GSE167186 <- read.csv('Count Data/GSE167186_counts.csv')
gse167186 <- getGEO('GSE167186', destdir = '.', getGPL = F)  
meta_167186 <- pData(phenoData(gse167186[[1]]))

## sifting out my columns of interest for conditions df
cond_167186 <- meta_167186 %>% 
  select(1,10,) %>%
  rename(muscle_health = characteristics_ch1) %>%
  mutate(muscle_health = gsub('group: ','', muscle_health)) %>% 
  mutate(muscle_health = gsub('Old Healthy','Old_Healthy', muscle_health)) %>%
  mutate(muscle_health = gsub('Young Healthy','Young_Healthy', muscle_health))
rownames(cond_167186) <- cond_167186$title
## duplicate count matrix
gene_count167186 <- count_GSE167186

## data manipulation to design deseqmtarix analysis 
rownames(gene_count167186) <- gene_count167186$Symbol
gene_count167186 <- subset(gene_count167186, select = -Symbol)
rownames(cond_167186) <- colnames(gene_count167186)

## create dds object 
dds_167186 <- DESeqDataSetFromMatrix(countData = gene_count167186, colData = cond_167186,
                        design = ~ muscle_health) 

## set factor levels
dds_167186$muscle_health <- factor(dds_167186$muscle_health, levels = c('Sarcopenia',
                                                                        'Old_Healthy', 'UNCLASSIFIED',  'Young_Healthy'))

## filering deseq data
keep_167186 <- rowSums(counts(dds_167186)) >= 10
dds_167186 <- dds_167186[keep,]

## run deseq
dds_167186 <- DESeq(dds_167186)
resultsNames(dds_167186)
## saving deseq as a result
res_167186 <- results(dds_167186)
res_167186
## variance standard transformation
vsd_167186 <- vst(dds_167186, blind = F)
vsd_167186

 


##--------------GSE151066----------
## DDS for GSE151066 
count_GSE151066 <- fread('Count Data/GSE151066_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count151066 <- as.data.frame(count_GSE151066, row.names = count_GSE151066$GeneID)
gene_count151066 <- gene_count151066[, -1]
gse151066 <- getGEO('GSE151066', destdir = '.', getGPL = F)
meta_151066 <- pData(phenoData(gse151066[[1]]))

## creating a coldata df
cond_151066 <- meta_151066 %>% 
  dplyr::select(1,10,11) %>%
  dplyr::rename(cohort = characteristics_ch1) %>%
  dplyr::rename(time_point = characteristics_ch1.1) %>%
  dplyr::mutate(cohort = gsub('cohort: ','', cohort)) %>%
  dplyr::mutate(time_point = gsub('time.point: ','', time_point))

##filtering out unneeded data from series matrix - everything that has 'post' in rows 
cond_151066 <- cond_151066[!cond_151066$time_point %in% c('Post','3hrPost') ,] 
 
## identify columns that are present as rownames in cond df - and filter count df by it
gene_count151066<- gene_count151066[, colnames(gene_count151066)  %in% rownames(cond_151066)]

## setting factor
cond_151066$cohort <- factor(cond_151066$cohort)

## dds object 
dds_151066 <- DESeqDataSetFromMatrix(countData = gene_count151066,
                                     colData = cond_151066, design = ~ cohort)

## filter dds by adhesome genes 
dds_151066 <- dds_151066[names(dds_151066) %in% ahhesome_genes$entrezgene_id,] 

## variance stabilization transformatio 
vsd_151066 <- vst(dds_151066, blind = F)
vsd_151066_df <- as.data.frame(assay(vsd_151066))
