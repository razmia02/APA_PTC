

############## SCRIPT FOR SCRNA-SEQ ANALYSIS OF GEO DATASET: GSE250521 ###############
getwd()

########## Loading the required libraries & pacakges #############################

library(Seurat)
library(dplyr)
library(SoupX)


############ STEP-1: LOAD THE DATASET ##########################

########## Define samples & paths ##################

sample_dir <- c("Data/scrna_data/GSM7980876_N-2", 
                "Data/scrna_data/GSM7980877_N-3", 
                "Data/scrna_data/GSM7980878_N-4", 
                "Data/scrna_data/GSM7980879_PTC-2", 
                "Data/scrna_data/GSM7980880_PTC-3", 
                "Data/scrna_data/GSM7980881_PTC-4")

sample_names <- c("Normal-1", "Normal-2", "Normal-3",
                  "PTC-1", "PTC-2", "PTC-3")


########## Load the data files, create seurat object & perform initial QC ###########


names(sample_dir) <- sample_names

list.files(sample_dir)

data <- Read10X(data.dir = sample_dir)

head(data)


dim(data) ## 60656 genes & 50172 cells


seurat <- CreateSeuratObject(counts = data, project = "PTC", 
                             min.cells = 3, min.features = 200)


####### Add the sample information (Normal or tumor) ##############


head(seurat$nFeature_RNA)

seurat$condition <- ifelse(grepl("Normal", seurat$orig.ident), "Normal", "PTC")

head(seurat$condition)



################ STEP-2: QUALITY CONTROL ################################


###### Calculate % of mito genes ###########

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT")

head(seurat@meta.data$percent.mt)

###### Plot the no of genes, transcripts and percent.mt in the raw data #########

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)


####### The plot shows some samples have less nfetaures than others ##########

####### Let's check the statistics per sample ############3

qc_summary <- seurat@meta.data  %>%
  group_by(orig.ident) %>%
  summarise(
    median_feature = median(nFeature_RNA), 
    median_counts = median(nCount_RNA), 
    median_mito = median(percent.mt),
    n_cells = n()
  )

qc_summary

##### The summary shows that PTC-3 has very high no of cells than other samples ####

