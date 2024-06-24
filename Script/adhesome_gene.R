
library(dplyr)
library(GEOquery)
library(tidyverse) 
library(DESeq2)
library(data.table)
library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

# install this to allow the ensmbl work

##-----adhesome gene db update-----
## choosing gene database to use 
ensembl<- useEnsembl(biomart = "genes") 
datasets <- listDatasets(ensembl) 
## selecting required dataset from database
ensemble_con <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = ensembl)
#building my updated adhesome gene index - here I mapped entrez and ensmbl gene ID to the adhesome df
attributes <- listAttributes(ensemble_con)
filters <- listFilters(ensemble_con)
adhesome_gene_map <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id', 'hgnc_id'), 
                  filters = "hgnc_id",
                  values = ahhesome_genes$HGNC_ID , mart = ensemble_con)

#merge adhesome gene to update with entrez and ensemble gene ID
ahhesome_genes<- ahhesome_genes %>%
  left_join(., adhesome_gene_map, by = c('HGNC_ID' = 'hgnc_id'))




##-----GSE126865------ 
## importing counts table and series matrix
count_GSE126865 <- fread("Count Data/GSE126865_raw_counts_GRCh38.p13_NCBI.tsv.gz")

count_GSE126865 <- as.data.frame(count_GSE126865, row.names = count_GSE126865$GeneID)
count_GSE126865 <- count_GSE126865[, -1]
gse126865 <- getGEO('GSE126865', destdir = '.', getGPL = F)
meta_126865 <- pData(phenoData(gse126865[[1]])) 

## creating coldata from series matrix metadata
cond_126865<- meta_126865 %>%
  dplyr::select(1,10) %>%
  dplyr::rename(group = characteristics_ch1) %>%
  mutate(group = gsub('group: before bed rest', 'before', group)) %>% 
  mutate(group = gsub('group: 10 days after bed rest', 'after', group))

all(colnames(count_GSE126865) == rownames(cond_126865)) 

## create dds object
dds_126865 <- DESeqDataSetFromMatrix(countData = count_GSE126865, 
                                     colData = cond_126865, design = ~ group)

## filter dds object by adhesome genes
dds_126865 <- dds_126865[names(dds_126865) %in% ahhesome_genes$entrezgene_id,]

## vsd
vsd_126865 <- vst(dds_126865, blind = F)
vsd_126865_df <- as.data.frame(assay(vsd_126865))
## save as csv file 
write.csv(vsd_126865_df, file = "GSE126865_variance.csv", row.names = T) 


  








##---GSE175495-----
## DDS for GSE175495
## importing counts and series matrix
count_GSE175495 <- fread('Count Data/GSE175495_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count175495 <- as.data.frame(count_GSE175495, row.names = count_GSE175495$GeneID)
gene_count175495 <- gene_count175495[, -1]

gse175495 <- getGEO('GSE175495', destdir = '.', getGPL = F)
meta_175495 <- pData(phenoData(gse175495[[1]]))

cond_175495 <- meta_175495 %>%
  dplyr::select(1,8,10) %>%
  dplyr::rename(source_name = source_name_ch1) %>%
  dplyr::rename(age = characteristics_ch1) %>%
  mutate(age = gsub('age: ','', age))

## filtering out unwanted samples from series coldata
cond_175495<- cond_175495[!cond_175495$source_name %in% c('Whole tissue (Adipose Tissue)_Young',
                                            'Whole tissue (Adipose Tissue)_Old'),]
## ensure colnames match rownames in coldata
gene_count175495 <- gene_count175495[, colnames(gene_count175495) %in% rownames(cond_175495)]
## incomplete match - code to match rownames with colnames
cond_175495 <- cond_175495[rownames(cond_175495) %in% colnames(gene_count175495), , drop = F]

all(colnames(gene_count175495) %in% rownames(cond_175495))
all(rownames(cond_175495) == colnames(gene_count175495))

## set factor 
cond_175495$age <- factor(cond_175495$age)

## DDS 
dds_175495 <- DESeqDataSetFromMatrix(countData = gene_count175495, 
                       colData = cond_175495, design = ~ age )
## filter by adhesome genes
dds_175495 <- dds_175495[names(dds_175495) %in% ahhesome_genes$entrezgene_id,]

## variance stabilization transformation
vsd_175495 <- vst(dds_175495, blind = F) 
vsd_175495_df <- as.data.frame(assay(vsd_175495))
 


##-----GSE113165------
## DDS for GSE113165
## import count and series matrix
count_GSE113165 <- fread('Count Data/GSE113165_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count113165 <- as.data.frame(count_GSE113165, row.names = count_GSE113165$GeneID)
gene_count113165 <- gene_count113165[, -1]

gse113165 <- getGEO('GSE113165', destdir = '.', getGPL = F)
meta_113165 <- pData(phenoData(gse113165[[1]]))
cond_113165 <- meta_113165 %>%
  select(1,11,12) %>%
  dplyr::rename(time = characteristics_ch1.1) %>% 
  dplyr::rename(age = characteristics_ch1.2) %>%
  mutate(age = gsub('age: ','', age)) %>%
  mutate(time = gsub('time: ', '', time)) %>%
  mutate(time = gsub('pre bed rest', 'pre_bed_rest', time)) %>%
  mutate(time = gsub('post bed rest', 'post_bed_rest', time))

## filtering and manipulating df
cond_113165<- cond_113165[!cond_113165$time %in% 'post_bed_rest',]
gene_count113165 <- gene_count113165[, colnames(gene_count113165) %in% rownames(cond_113165)]

all(colnames(gene_count113165) == rownames(cond_113165))
all(rownames(cond_113165) %in% colnames(gene_count113165))

##factor
cond_113165$age <- factor(cond_113165$age)
## DDS object
dds_113165 <- DESeqDataSetFromMatrix(countData = gene_count113165, 
                       colData = cond_113165, design = ~ age) 
## filter by adhesome genes
dds_113165 <- dds_113165[names(dds_113165) %in% ahhesome_genes$entrezgene_id,]

## variance stabilization trans
vsd_113165 <- vst(dds_113165, blind = F) 
vsd_113165_df <- as.data.frame(assay(vsd_113165)) 

 
 


##---GSE164471-----
## DDS for GSE164471
## import counts and series matrix
count_GSE164471 <- fread('Count Data/GSE164471_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count164471 <- as.data.frame(count_GSE164471, row.names = count_GSE164471$GeneID)
gene_count164471 <- gene_count164471[, -1]
gse164471 <- getGEO('GSE164471', destdir = '.', getGPL = F) 
meta_164471 <- pData(phenoData(gse164471[[1]]))

## coldata design
cond_164471 <- meta_164471 %>%
  select(1, 11,12) %>%
  dplyr::rename(age = characteristics_ch1.1)%>%
  dplyr::rename(gender = characteristics_ch1.2) %>%
  mutate(age = gsub('age: ', '', age)) %>%
  mutate(gender = gsub('gender: ', '', gender)) 

all(colnames(gene_count164471) == rownames(cond_164471))
all(rownames(cond_164471) %in% colnames(gene_count164471))
## factor set
cond_164471$age <- factor(cond_164471$age)
##DDS
dds_164471 <- DESeqDataSetFromMatrix(countData = gene_count164471, colData = cond_164471,
                       design = ~ age)
## filter by adhesome genes 
dds_164471 <- dds_164471[names(dds_164471) %in% ahhesome_genes$entrezgene_id,]

## variance stab trans
vsd_164471 <- vst(dds_164471, blind = F)
vsd_164471_df <- as.data.frame(assay(vsd_164471))





##----GSE165630----
## DDS for GSE165630
count_GSE165630 <- fread('Count Data/GSE165630_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count165630 <- as.data.frame(count_GSE165630, row.names = count_GSE165630$GeneID)
gene_count165630 <- gene_count165630[, -1]
gse165630 <- getGEO('GSE165630', destdir = '.', getGPL = F)
meta_165630 <- pData(phenoData(gse165630[[1]]))
cond_165630 <- meta_165630 %>%
  select(1,11) %>%
  dplyr::rename(age = characteristics_ch1.1) %>%
  mutate(title = gsub('aged.*', '', title)) %>%
  mutate(title = gsub('muscle of *', '', title)) %>%
  mutate(title = gsub('well ', '', title)) %>%
  mutate(title = gsub('resistance trained', 'resistance_trained', title)) %>%
  mutate(title = gsub('endurance trained ', 'endurance_trained', title)) %>%
  mutate(age = gsub('age \\(years\\): ', '', age))

all(rownames(cond_165630) == colnames(gene_count165630))
all(colnames(gene_count165630) %in% rownames(cond_165630))

## factor set
cond_165630$title <- factor(cond_165630$title)
## DDS 
dds_165630 <- DESeqDataSetFromMatrix(countData = gene_count165630, 
                       colData = cond_165630, design = ~1 + title) 

## filter by adhesome genes
dds_165630 <- dds_165630[names(dds_165630) %in% ahhesome_genes$entrezgene_id,]

## variance standard trans
vsd_165630 <- vst(dds_165630, blind = F)
vsd_165630_df <- as.data.frame(assay(vsd_165630))
 


##----GSE120642----
## DDS for GSE120642
count_GSE120642 <- fread('Count Data/GSE120642_raw_counts_GRCh38.p13_NCBI.tsv.gz')
gene_count120642 <- as.data.frame(count_GSE120642, row.names = count_GSE120642$GeneID)
gene_count120642 <- gene_count120642[, -1]

gse120642 <- getGEO('GSE120642', destdir = '.', getGPL = F)
meta_120642 <- pData(phenoData(gse120642[[1]]))
cond_120642 <- meta_120642 %>%
  select(1,11) %>%
  dplyr::rename(diagnosis = characteristics_ch1.1) %>%
  mutate(diagnosis = gsub('diagnosis: ','', diagnosis)) %>%
  mutate(diagnosis = gsub('Healthy Adult','Healthy_Adult', diagnosis)) %>%
  mutate(diagnosis = gsub('Intermittent claudicant','Intermittent_claudicant', diagnosis)) %>%
  mutate(diagnosis = gsub('Critical limb ischemia','Critical_limb_ischemia', diagnosis))

all(colnames(gene_count120642) == rownames(cond_120642))
all(rownames(cond_120642) %in% colnames(gene_count120642))
## set factor
cond_120642$diagnosis <- factor(cond_120642$diagnosis)
## DDS
dds_120642 <- DESeqDataSetFromMatrix(countData = gene_count120642, 
                       colData = cond_120642, design = ~ 1 + diagnosis)

## filter by adhesome genes
dds_120642 <- dds_120642[names(dds_120642) %in% ahhesome_genes$entrezgene_id,]

## variance 
vsd_120642 <- vst(dds_120642, blind = F)
vsd_120642_df <- as.data.frame(assay(vsd_120642))
 
