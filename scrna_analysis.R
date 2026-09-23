

############## SCRIPT FOR SCRNA-SEQ ANALYSIS OF GEO DATASET: GSE250521 ###############
getwd()

########## Loading the required libraries & pacakges #############################

library(Seurat)
library(dplyr)
library(SoupX)
library(remotes)
library(knitr)
library(rmarkdown)
library(scQCenrich)
library(mclust)
library(scDblFinder)
library(singlecellexpreiment)


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



######### Lets try to run scQCenrich to understand per sample metrics ##########

list_panglao_tissues() ### Check the list of tissues available in db

####### Subset the PTC-3 sample ###########

# table(seurat$orig.ident)
# 
# ptc_3 <- subset(seurat, subset = orig.ident == "PTC-3")
# 
# 
# qc_ptc_3 <- run_qc_pipeline(
#   obj         = ptc_3,
#   species     = "human",
#   tissue      = c("Thyroid"),
#   method      = "gmm",
#   qc_strength = "auto",
#   report_file = "Results/scQCenrich_ptc-3.html")
# 
# head(qc_ptc_3$status_df)


############### Loop through the process ###################

######## 1. Subset samples from seurat object ##################

######## 2. Run scQCenrich on each sample ############

######## 3. Select the cells to keep (as decided by scQCenrich) #######

sample_list <- unique(seurat$orig.ident)

qc_status_list <- list()

for (s in sample_list) {
  
  sub_obj <- subset(seurat, subset = orig.ident == s)
  
  res <- run_qc_pipeline(
    obj         = sub_obj,
    species     = "human",
    tissue      = c("Thyroid"),
    method      = "gmm",
    qc_strength = "auto",
    report_file = paste0("Results/scQCenrich_", s, ".html")
  )
  
  qc_status_list[[s]] <- res$status_df
}


########## Remove the Doublets from each sample #####################

####### 1. Select the cells to keep from scQCenrich results ###########

####### 2. Run scDblFinder to identify doublets ##########

####### 3. Keep the singlets ##################


