BiocManager::install("SingleCellExperiment")
install.packages("SingleR")
install.packages("SingleR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")

install.packages("rmarkdown", dependencies = TRUE)


.libPaths("C:\\Program Files\\R\\R-4.1.2\\library")

library(Seurat)
library(patchwork)
library(dplyr)
library(scran)
library(srt)
library(ggplot2)
library(cowplot)
library(celldex)
library(SingleR)
library(scmap)
library(SingleCellExperiment)

# Load the wt dataset 5 weeks
wt.data <- Read10X(data.dir = "E:/AMIN HDD/STUDY MATERIALS/UNIVERSITY OF TSUKUBA/R Coding/Fibrillin_1/5 Weeks Old Homozygous and WT mice/5 Weeks old WT/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
wt <- CreateSeuratObject(counts = wt.data, project = "5 weeks WT", min.cells = 3, min.features = 200)
wt

# Load the homozygous data set 5 weeks
homozygous.data <- Read10X(data.dir = "E:/AMIN HDD/STUDY MATERIALS/UNIVERSITY OF TSUKUBA/R Coding/Fibrillin_1/5 Weeks Old Homozygous and WT mice/5 weeks old homozygous/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
homozygous <- CreateSeuratObject(counts = homozygous.data, project = "5 Weeks Homozygous", min.cells = 3, min.features = 200)
homozygous

# Load the homozygous data set 3 weeks
homozygous_3.data <- Read10X(data.dir = "E:/AMIN HDD/STUDY MATERIALS/UNIVERSITY OF TSUKUBA/R Coding/Fibrillin_1/3 Weeks Old Homozygous and WT data/3 weeks homozygous")
# Initialize the Seurat object with the raw (non-normalized data).
homozygous_3 <- CreateSeuratObject(counts = homozygous_3.data, project = "3 Weeks Homozygous", min.cells = 3, min.features = 200)
homozygous_3


# Create Seurat Object 
wt   <- CreateSeuratObject(wt.data, project = "wt")
homozygous   <- CreateSeuratObject(homozygous.data, project = "homozygous_5")
homozygous_3   <- CreateSeuratObject(homozygous_3.data, project = "homozygous_3")

# Remove matrices to save memory
rm(wt.data)
rm(homozygous.data)
rm(homozygous_3.data)

# Calculate the fractions of mitochondrial genes and ribosomal proteins, and do quick-and-dirty filtering of the datasets
wt[["percent.mt"]]  <- PercentageFeatureSet(wt, pattern = "^mt-")
wt[["percent.rbp"]] <- PercentageFeatureSet(wt, pattern = "^RP[SL]")

homozygous[["percent.mt"]]  <- PercentageFeatureSet(homozygous, pattern = "^mt-")
homozygous[["percent.rbp"]] <- PercentageFeatureSet(homozygous, pattern = "^RP[SL]")

homozygous_3[["percent.mt"]]  <- PercentageFeatureSet(homozygous_3, pattern = "^mt-")
homozygous_3[["percent.rbp"]] <- PercentageFeatureSet(homozygous_3, pattern = "^RP[SL]")


VlnPlot(wt, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

VlnPlot(homozygous, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)

VlnPlot(homozygous_3, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)



table(rownames(wt) %in% rownames(homozygous) %in% rownames(homozygous_3)) 

table(rownames(wt) %in% rownames(homozygous)) #%in% rownames(homozygous_3)) 


# Quick filtering of the datasets removes dying cells and putative doublets
wt <- subset(wt, subset = nFeature_RNA >= 200 & percent.mt <= 10)
homozygous <- subset(homozygous, subset = nFeature_RNA >= 200 & percent.mt <= 10)
homozygous_3 <- subset(homozygous_3, subset = nFeature_RNA >= 200 & percent.mt <= 10)


wt
homozygous
homozygous_3

combined <- list()

combined [["wt"]] <- wt
combined [["homozygous"]] <- homozygous
combined [["homozygous_3"]] <- homozygous_3

for (i in 1:length(combined)) {
  combined[[i]] <- NormalizeData(combined[[i]], verbose = F)
  combined[[i]] <- FindVariableFeatures(combined[[i]], selection.method = "vst",nfeatures = 2000, verbose = F)
}

combined_anchors    <- FindIntegrationAnchors(object.list = combined, dims = 1:20)

combined_seurat     <- IntegrateData(anchorset = combined_anchors, dims = 1:20)

rm(combined)
rm(combined_anchors)

DefaultAssay(combined_seurat) <- "RNA"

combined_seurat <- NormalizeData(combined_seurat, verbose = F)

combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 2000, verbose = F)

combined_seurat <- ScaleData(combined_seurat, verbose = F)

combined_seurat <- RunPCA(combined_seurat, npcs = 50, verbose = F)

combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:20, verbose = F)

ElbowPlot(combined_seurat, ndims = 20, reduction = "pca")

DimPlot(combined_seurat, reduction = "umap") + plot_annotation(title = "before integration of sham and aaa datasets")

DefaultAssay(combined_seurat) <- "integrated"

combined_seurat <- ScaleData(combined_seurat, verbose = F)

combined_seurat <- RunPCA(combined_seurat, npcs = 50, verbose = F)

combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:20, verbose = F)

DimPlot(combined_seurat, reduction = "umap") + plot_annotation(title = "integration of WT and Homo")

DimPlot(combined_seurat, reduction = "umap", split.by = "orig.ident")

combined_seurat <- FindNeighbors(combined_seurat, dims = 1:20, k.param = 70, verbose = F)

combined_seurat <- FindClusters(combined_seurat, verbose = F)


DimPlot(combined_seurat,label = T) 

count_table <- table(combined_seurat@meta.data$seurat_clusters, combined_seurat@meta.data$orig.ident)
count_table
# Determine the "nearest neighbors" of each cell
combined_seurat <- Seurat::FindNeighbors(combined_seurat, dims = 1:15)
# Cluster the cells
combined_seurat <- Seurat::FindClusters(combined_seurat, resolution = 0.5)
# visualize our data on a UMAP
Seurat::DimPlot(combined_seurat, reduction = "umap", label = T)

# Now let's extract the top marker genes, and see which ones correspond with each cluster. This can be done using the FindAllMarkers function within Seurat
markers_seur <- Seurat::FindAllMarkers(combined_seurat, only.pos = TRUE)
# to show the marker genes
markers_seur
# Retrieve the top 20 marker genes per cluster
top20 <- markers_seur %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_seur), value = TRUE)),
                   n = 20)

top20


# check the marker genes in each clusters

FeaturePlot(combined_seurat, features = c("Pecam1", "Cdh5", "Cldn5")) 

FeaturePlot(combined_seurat, features = c("Cd3d", "Cd3g","Cd28")) 

FeaturePlot(combined_seurat, features = c("Cd79a", "Ly6d","Cd79b")) 

FeaturePlot(combined_seurat, features = c("Klrd1", "Flt3","H2-Ab1",)) 
FeaturePlot(combined_seurat, features = c("Cd3g", "Gzma"))
FeaturePlot(combined_seurat, features = c("S100a8", "S100a9", "Cxcr4")) 
FeaturePlot(combined_seurat, features = c("Myh11", "Tagln", "Acta2")) 

FeaturePlot(combined_seurat, features = c("Dcn", "Lum", "Col1a1"))

FeaturePlot(combined_seurat, features = c("Lyz2", "Pf4", "Mrc1", "Il1b")) 

FeaturePlot(combined_seurat, features = c("Cd68", "Cdca3", "C1qb")) 

FeaturePlot(combined_seurat, features = c("H2-Ab1", "Stmn1", "Mki67"))
FeaturePlot(combined_seurat, features = c("Neat1", "Cd79a"))


plots <-VlnPlot(combined_seurat, features = c("Neat1"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)


# Put new cluster ID after identification
new.cluster.ids <- c("SMC-1", "SMC-2", "Macro-1", "SMC-3", "SMC-4", "Fibro-1", "Marco-3","Neutro", "Macro-4", "NK & T Cells", "SMC-5", "SMC-6", "ECs")
names(new.cluster.ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, new.cluster.ids)
Seurat::DimPlot(combined_seurat, reduction = "umap", label = T)


# count number of cells in individual clusters 
count_table <- table(combined_seurat@meta.data$seurat_clusters, combined_seurat@meta.data$orig.ident)
count_table

# Now let's extract the top marker genes, and see which ones correspond with each cluster. This can be done using the FindAllMarkers function within Seurat
markers_seur <- Seurat::FindAllMarkers(combined_seurat, only.pos = TRUE)
# to show the marker genes
markers_seur
# Retrieve the top 5 marker genes per cluster
top5 <- markers_seur %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(markers_seur), value = TRUE)),
                   n = 5)

top5

# Plot the cells using dot plot after Cell-type annotation
Seurat::DotPlot(combined_seurat, features = unique(top5$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) 
# Plot the cells using Heatmap after Cell-type annotation
Seurat::DoHeatmap(combined_seurat, features = unique(top5$gene)) +
  Seurat::NoLegend() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

# Cluster 0 Expression
plots <-VlnPlot(combined_seurat, features = c("Tagln", "Acta2", "Tpm2", "Myl9", "Myh11"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Tagln", "Acta2"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Tpm2", "Myl9"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Myh11"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 1 expression
plots <-VlnPlot(combined_seurat, features = c("Ccl12", "Hspa1b", "Pf4", "Hspa1a", "Mgp"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Ccl12", "Hspa1b"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Pf4", "Hspa1a"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Mgp"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 2 expression
plots <-VlnPlot(combined_seurat, features = c("C1qa", "C1qc", "C1qb", "Ccl7", "Ccl2"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("C1qa", "C1qc"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("C1qb", "Ccl7"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Ccl2"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 3 expression
plots <-VlnPlot(combined_seurat, features = c("Dcn", "Fbln1", "Col1a1", "Serpinf1", "Col1a2"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Dcn", "Fbln1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Col1a1", "Serpinf1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Col1a2"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 4 expression
plots <-VlnPlot(combined_seurat, features = c("Col1a2", "S100a8", "Retnlg", "Srgn", "Slpi"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Col1a2", "S100a8"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Retnlg", "Srgn"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Slpi"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))


# Cluster 5 expression
plots <-VlnPlot(combined_seurat, features = c("Lyz2", "Plac8", "Chil3", "Lgals3", "Cybb"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Lyz2", "Plac8"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Chil3", "Lgals3"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Cybb"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 6 expression
plots <-VlnPlot(combined_seurat, features = c("Cd74", "H2-Eb1", "H2-Aa", "H2-Ab1", "Cytip"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Cd74", "H2-Eb1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("H2-Aa", "H2-Ab1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Cytip"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 7 expression
plots <-VlnPlot(combined_seurat, features = c("Trdc", "Cd3g", "Emb", "Ctla2a", "Il5"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Trdc", "Cd3g"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Emb", "Ctla2a"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Il5"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))


# Cluster 8 expression
plots <-VlnPlot(combined_seurat, features = c("Cytl1", "Pecam1", "Edn1", "Vcam1", "Ptprb"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Cytl1", "Pecam1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Edn1", "Vcam1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Ptprb"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 9 expression
plots <-VlnPlot(combined_seurat, features = c("Igfbp5", "Pi16", "Nid1", "C3", "Smpd3"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Igfbp5", "Pi16"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Nid1", "C3"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("C3", "Smpd3"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

# Cluster 10 expression
plots <-VlnPlot(combined_seurat, features = c("Rgs5", "Eln", "Smoc1", "Htra1", "Sod3"), split.by = "orig.ident", pt.size = 0)
wrap_plots(plots = plots, ncol = 1)

FeaturePlot(combined_seurat, features = c("Rgs5", "Eln"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Smoc1", "Htra1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(combined_seurat, features = c("Sod3"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))





query_seur <- cerebroApp::getMarkerGenes(combined_seurat,
                                         groups = c('seurat_clusters'),
                                         assay = "RNA",
                                         organism = "mm", samples2compare="all", 
                                         annotate=TRUE, 
                                         chip=NULL,
                                         score.cutoff=1)

query_seur <- cerebroApp::getEnrichedPathways(combined_seurat,
                                              databases = c("GO_Biological_Process_2018",
                                                            "GO_Cellular_Component_2018",
                                                            "GO_Molecular_Function_2018",
                                                            "KEGG_2016",
                                                            "WikiPathways_2016",
                                                            "Reactome_2016",
                                                            "Panther_2016",
                                                            "Human_Gene_Atlas",
                                                            "Mouse_Gene_Atlas"),
                                              adj_p_cutoff = 0.05,
                                              max_terms = 100,
                                              URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich")



