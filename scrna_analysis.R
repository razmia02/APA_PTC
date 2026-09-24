

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
library(SingleCellExperiment)
library(ggplot2)


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

######## Combine all sample's qc into one table ###########

qc_status_merge <- do.call(rbind, unname(qc_status_list))


colnames(qc_status_list[["PTC-3"]]) #### Check the colnames 

##### qc_status contains data to keep/remove/boderline cells ########

######### ADD QC_STATUS OF CELLS TO SEURAT OBJECT ################

seurat$qc_status <- qc_status_merge[colnames(seurat), "qc_status"]

table(seurat$qc_status) ##### sanity check if the qc status is merged


########## Remove the Doublets from each sample #####################

####### 1. Select the cells to keep from scQCenrich results ###########

####### 2. Run scDblFinder to identify doublets ##########

####### 3. Keep the singlets ##################


sample_names <- unique(seurat$orig.ident)

filtered_list <- list()

dbl_status_list <- list()

for (s in sample_names) {
  sub_obj <- subset(seurat, subset = orig.ident == s & qc_status == "keep")
  
  ##### Convert the seurat object into singlecellexperiment
  
  sce <- as.SingleCellExperiment(sub_obj)
  
  sce <- scDblFinder(sce)
  
  #### add the scdblfinder results back to seurat object 
  
  sub_obj$scDblFinder.class <- sce$scDblFinder.class
  
  #### sanity check 
  
  print(table(sub_obj$scDblFinder.class))
  
  #### save the proportions of doublet/singlet per sample
  
  dbl_status_list[[s]] <- data.frame(
    
    orig.ident = sub_obj$orig.ident,
    
    scDblFinder.class = sub_obj$scDblFinder.class
  )
  
  #### subset to singlets 
  
  sub_obj_sin <- subset(sub_obj, subset = scDblFinder.class == "singlet")
  
  filtered_list[[s]] <- sub_obj_sin
}

dbl_status_df <- do.call(rbind, unname(dbl_status_list))

################# Merge the filtered list ##########################

seurat_filt <- merge(filtered_list[[1]], y = filtered_list[-1], 
                     project = "APA_thyroid")

seurat_filt <- JoinLayers(seurat_filt)

table(seurat_filt$orig.ident) #### sanity check

table(seurat_filt$scDblFinder.class) #### sanity check


######### Visualize the doublets per sample #################

ggplot(dbl_status_df, aes(x = orig.ident, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(y = "Proportion", title = "Doublet Proportion by Sample")


####### Summary of doublet/singlet proportions per sample ##########

dbl_summary_wide <- dbl_status_df %>%
  group_by(orig.ident, scDblFinder.class) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    id_cols = orig.ident,
    names_from = scDblFinder.class,
    values_from = c(n, percent)
  )

dbl_summary_wide

##### PTC-3 have the highest proportion of doublets #############


########## Normalization ################3

