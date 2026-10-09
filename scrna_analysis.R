

############## SCRIPT FOR SCRNA-SEQ ANALYSIS OF GEO DATASET: GSE250521 ###############
getwd()

########## Loading the required libraries & pacakges #############################

library(Seurat)
library(dplyr)
library(scQCenrich)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(harmony)
library(copykat)
library(devtools)

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


########## FINAL SEURAT OBJECT: seurat_filt ####################

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


########## Normalization ################

seurat_filt <- NormalizeData(seurat_filt)


############### Find Variable Features #####################

seurat_filt <- FindVariableFeatures(seurat_filt, nfeatures = 3000)


top_features <- head(VariableFeatures(seurat_filt), 20)

plot1 <- VariableFeaturePlot(seurat_filt)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2


############### Scale Data #########################3

seurat_filt <- ScaleData(seurat_filt)


############## Linear Dimensionality Reduction (PCA) #####################

seurat_filt <- RunPCA(seurat_filt, npcs = 50)

ElbowPlot(seurat_filt, ndims = ncol(Embeddings(seurat_filt, "pca")))

###### Elbow plot shows first 30-35 PCs showing the variation ###########

##### Lets select first 15 PCs for downstream analysis ############

############ Batch Correction using Harmony ################

###### Run UMAP before & after to compare batch correction results ############

seurat_filt <- RunUMAP(seurat_filt, reduction = "pca", 
                  dims = 1:15, 
                  reduction.name = "umap.unintegrated")

plot1 <- DimPlot(seurat_filt, reduction = "umap.unintegrated", 
        group.by = "condition") + ggtitle("Before Harmony")


######### Run Harmony #########

seurat_filt <- RunHarmony(seurat_filt, 
                     group.by.vars = "condition", 
                     dims.use = 1:15, max_iter = 50)


########### Non-Linear Dimensionality Reduction (UMAP) ################

seurat_filt <- RunUMAP(seurat_filt, reduction = "harmony", 
                  dims = 1:15, reduction.name = "umap")


plot2 <- DimPlot(seurat_filt, reduction = "umap", 
                 group.by = "condition") + ggtitle("After Harmony")

plot2


plot1 + plot2 #### See if the batch correction was applied correctly


########### Find Neighbours #################

seurat_filt <- FindNeighbors(seurat_filt, dims = 1:15)

############ Find Clusters ##################

seurat_filt <- FindClusters(seurat_filt, resolution = 0.3) #### 17 clusters


DimPlot(seurat_filt, reduction = "umap", label = TRUE)


########### Find cluster markers ################

cl_markers <- FindAllMarkers(seurat_filt, only.pos = TRUE, min.pct = 0.25, 
                             logfc.threshold = log(1.2))

View(cl_markers)

cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


############ Annotate clusters #################


thyroid_markers <- list(
  Thyrocyte    = c("TG", "TPO", "PAX8", "TSHR", "KRT8", "KRT18", "NKX2-1"),
  B_cell       = c("CD19", "CD79A", "MS4A1", "CD20","CD79B", "IGKC", "CD74", "CD22", "CD27"),
  Myeloid      = c("LYZ", "FCER1G", "TYROBP", "CD68", "CD14", "C1QA", "CD163", "CD33"),
  T_cell       = c("CD3", "CD3D", "CD3E", "CD4", "CD8", "CD2", "CD8A", "CD69", "CD28"),
  Endothelial  = c("PECAM1", "VWF", "CLDN5", "CDH5", "CD34", "KDR", "FLT1"),
  Fibroblast   = c("FAP", "DCN", "ACTA2", "COL1A1", "COL1A2", "DCN", "PDGFRA"),
  SmoothMuscle = c("ACTA2", "RGS5", "MYH11", "TAGLN", "CNN1", "MYL9"),
  Proliferating_cell = c("MKI67", "TOP2A", "TYMS", "DUT")
)


VlnPlot(seurat_filt, features = c("CD11C", "CD1C", "CD83", "CD86", "FCER1A", "CD80", "CD74", "FLT3"), 
        pt.size = 0)


######## Cluster: 2,3,5,8 = thyrocytes ##############

######## Cluster: 12 = B_cell ###############

######## Cluster: 4,10 = Macrophage #########

######## Cluster: 0,11 = CD8+ T_cell ###########

######## Cluster: 1 = CD4+ T_cell ###########

######## Cluster: 9,16 = Endothelial_cell ############

######## Cluster: 13 = Fibroblast ###########

######## Cluster: 6 = Smooth_Muscle ##########

######## Cluster: 15 = Proliferating_Myeloid #########

######## Cluster: 7 = Unknown #########

cluster_7 <- FindMarkers(seurat_filt, ident.1 = "7")

View(cluster_7)


############# ANNOTATE THE CLUSTERS ###################


new_ident <- setNames(c("0" = "CD8+T_cell", 
                        "1" = "CD4+T_cell", 
                        "2" = "Thyrocyte", 
                        "3" = "Thyrocyte",
                        "4" = "Macrophage", 
                        "5" = "Thyrocyte", 
                        "6" = "Smooth_Muscle", 
                        "7" = "Unknown_1",
                       "8" = "Thyrocyte", 
                       "9" = "Endothelial_cell", 
                       "10" = "Macrophage", 
                       "11" = "CD8+T_cell",
                       "12" = "B_cell", 
                       "13" = "Fibroblast", 
                       "14" = "Unknown_2", 
                       "15" = "Proliferating_Myeloid",
                      "16"  = "Endothelial_cell"),
                      levels(seurat_filt))

seurat_filt <- RenameIdents(seurat_filt, new_ident)

DimPlot(seurat_filt, reduction = "umap", label = TRUE)

# 0 = CD8+ T_cell
# 1 = CD4+ T_cell
# 2 = Thyrocyte
# 3 = Thyrocyte
# 4 = Macrophage
# 5 = Thyrocyte
# 6 = Smooth_Muscle
# 7 = Unknown
# 8 = Thyrocyte
# 9 = Endothelial_cell
# 10 = Macrophage
# 11 = CD8+ T_cell
# 12 = B_cell
# 13 = Fibroblast
# 14 = Unknown
# 15 = Proliferating_Myeloid
# 16 = Endothelial_cell

saveRDS(seurat_filt, file = "seurat.rds")

seurat_filt <- readRDS("seurat.rds")


############ Subset the seurat object to one sample only ##################

cells <- colnames(seurat_filt)[seurat_filt$orig.ident == "Normal-1"]

length(cells)

mat <- LayerData(seurat_filt, assay = "RNA", layer = "data")[, cells]

rownames(mat)

################# Test run copykat ########################

copykat.test <- copykat(rawmat = mat, 
                        id.type="S", 
                        ngene.chr=5, 
                        win.size=25, 
                        KS.cut=0.1, 
                        sam.name="test", 
                        distance="euclidean", 
                        norm.cell.names="",
                        output.seg="FLASE", 
                        plot.genes="TRUE", 
                        genome="hg20")

################# ADDITIONAL CODE #######################

######### CHCEK HOW MANY PCS TO USE FOR DOWNSTREAM #################

# pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
# cumu <- cumsum(pct)
# 
# co1 <- which(cumu > 90 & pct < 5)[1]  # cumulative variance >90% AND this PC contributes <5% alone
# co2 <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
# # last point where consecutive PCs still differ by >0.1%
# 
# pcs_selected <- min(co1, co2)
# pcs_selected
