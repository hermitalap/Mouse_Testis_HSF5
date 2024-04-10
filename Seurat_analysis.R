## Loading dependent packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(scCustomize)
library(gridExtra)
library(grid)


## Read data
setwd("path/to/data/")
Hsf5.ko.data <- readRDS("Hsf5.ko.merge.raw.Rds")
Hsf5.wt.data <-readRDS("Hsf5.wt.merge.raw.Rds")

# Preprocessing
## Merge samples and add identity meta information
Hsf5.merge <- c(CreateSeuratObject(counts = Hsf5.ko.data, min.cells = 3, min.features = 200, assay = "RNA"),CreateSeuratObject(counts = Hsf5.wt.data, min.cells = 3, min.features = 200, assay = "RNA"))
Hsf5.merge = merge(Hsf5.merge[[1]], y = Hsf5.merge[-1], add.cell.ids = c("Hsf5.KO", "WT"))
orig.ident = rep("Hsf5.KO",ncol(Hsf5.ko.data))
orig.ident = append(orig.ident,rep("Hsf5.WT",ncol(Hsf5.wt.data)))
Hsf5.merge[["orig.ident"]] = orig.ident

sample = rep("Knockout",ncol(Hsf5.ko.data))
sample = append(sample,rep("Wildtype",ncol(Hsf5.wt.data)))
Hsf5.merge[["sample"]] = sample
Hsf5.merge[["percent.mt"]] <- PercentageFeatureSet(Hsf5.merge, pattern = "^mt-")
png("./Hsf5.merge_nFeature_and_percent.mt_Vln.png",res = 600, width=10, height=6, units="in")
VlnPlot(Hsf5.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()
Hsf5.merge <- subset(Hsf5.merge, subset = nFeature_RNA > 200 & nFeature_RNA < 8500 & percent.mt < 20)

# Routine analysis
Hsf5.merge <- NormalizeData(Hsf5.merge, normalization.method = "LogNormalize", scale.factor = 10000)
Hsf5.merge <- FindVariableFeatures(Hsf5.merge, selection.method = "vst", nfeatures = 2000)
Hsf5.merge <- ScaleData(Hsf5.merge)
Hsf5.merge <- RunPCA(Hsf5.merge, features = VariableFeatures(object = Hsf5.merge))
Hsf5.merge <- FindNeighbors(Hsf5.merge, dims = 1:15, reduction = "pca")
Hsf5.merge <- FindClusters(Hsf5.merge, resolution = 0.5, cluster.name = "unintegrated_clusters")
Hsf5.merge <- RunUMAP(Hsf5.merge, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(Hsf5.merge, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

## Harmony method removes batch effects
Hsf5.merge <- IntegrateLayers(object = Hsf5.merge, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",verbose = FALSE)
Hsf5.merge <- FindNeighbors(Hsf5.merge,  reduction = "harmony", dims = 1:15)
Hsf5.merge <- RunUMAP(Hsf5.merge, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
Hsf5.merge <- FindClusters(Hsf5.merge, resolution = 0.1,  cluster.name = "harmony_clusters")
DimPlot(Hsf5.merge,reduction = "umap.harmony",group.by = c("orig.ident", "harmony_clusters"))
png("./Raw UMAP analysis.png",res = 600, width=15, height=7, units="in")
DimPlot(Hsf5.merge,reduction = "umap.harmony",group.by = c("orig.ident", "harmony_clusters"))
dev.off()

# Analysis of whole testicular cells
Hsf5.merge.ident <- Hsf5.merge # make a copy

# all.markers <- FindAllMarkers(object = Hsf5.merge.ident, node = 4)
# write.csv(all.markers,"./data/Hsf5.merge.ident.markers_all.csv")
# FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Sox9"), order = T)
# Refer to known literature to determine the identity of each cell group

DimPlot(Hsf5.merge.ident, reduction = "umap.harmony",label = TRUE, group.by = c("orig.ident", "harmony_clusters"))
new.cluster.ids <- c('Leydig cells', 'Interstitial progenitors', 'SCytes', 'SPG', 'STids', 'SCytes', 'STids', 'Sertoli cells', 'STids', 'Sertoli cells', 'SPG', 'SCytes', 'SCytes', 'Macrophages', 'Endothelial cells', 'Perivascular cells', 'Peritubular myoid cells', 'Macrophages', 'STids', 'Interstitial progenitors')
names(new.cluster.ids) <- levels(Hsf5.merge.ident)
Hsf5.merge.ident <- RenameIdents(Hsf5.merge.ident, new.cluster.ids)
DimPlot(Hsf5.merge.ident, reduction = "umap.harmony",label = TRUE)   

write_rds(Hsf5.merge.ident, "../data/Whole_testis.rds")


# Extract germ cells
# Hsf5.merge.germ.cell <- Hsf5.merge.ident[,!Hsf5.merge.ident@active.ident %in% c('SPG', 'SCytes', 'STids')]
Hsf5.merge.germ.cell <- Hsf5.merge[,Hsf5.merge@meta.data$harmony_clusters %in% c(2,3,4,5,6,8,10,11,12,18)]
DimPlot(Hsf5.merge.germ.cell,reduction = "umap.harmony",label = TRUE, group.by = c("orig.ident", "harmony_clusters"))

Hsf5.merge.germ.cell <- FindNeighbors(Hsf5.merge.germ.cell,  reduction = "harmony", dims = 1:19)
Hsf5.merge.germ.cell <- RunUMAP(Hsf5.merge.germ.cell, reduction = "harmony", dims = 1:19, reduction.name = "umap.harmony", return.model=TRUE)
Hsf5.merge.germ.cell <- FindClusters(Hsf5.merge.germ.cell, resolution = 0.5,  cluster.name = "harmony_clusters")

DimPlot(Hsf5.merge.germ.cell,reduction = "umap.harmony",label = TRUE, group.by = c("orig.ident", "harmony_clusters"))
png("./Germ_cells_unidentified.png",res = 600, width=15, height=7, units="in")
DimPlot(Hsf5.merge.germ.cell,reduction = "umap.harmony",label = TRUE, group.by = c("orig.ident", "harmony_clusters"))
dev.off()

FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("nFeature_RNA","percent.mt"))
png("./Germ_cells_unidentified_features_and_percent_mt.png",res = 600, width=15, height=7, units="in")
VlnPlot(Hsf5.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

# Refer to known literature to determine the identity of each cell group
png("./Spermatogenesis_DEGs/1.FeaturePlot_spg.png",res = 600, width=15, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Zbtb16","Sohlh1","Kit","Uchl1","Dmrt1","Dazl","Stra8","Scml2","Dmc1"))
dev.off()

png("./Spermatogenesis_DEGs/1.FeaturePlot_Undiff.png",res = 600, width=10, height=10, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Egr4","Sox4","Sohlh1"))
dev.off()

png("./Spermatogenesis_DEGs/1.FeaturePlot_Diff1.png",res = 600, width=10, height=10, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Uchl1","Stra8","Dmrt1"))
dev.off()

png("./Spermatogenesis_DEGs/1.FeaturePlot_Diff2.png",res = 600, width=10, height=10, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Ccnd1","Plk1","Esx1"))
dev.off()

png("./Spermatogenesis_DEGs/2.FeaturePlot_L.png",res = 600, width=20, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Prss50","Smc3","Tex15","Smc1b","Syce2","Dmc1","Prdm9","Fbxw2","Brca2","Rec8","Arl6ip1","Nol8"))
dev.off()

png("./Spermatogenesis_DEGs/3.FeaturePlot_Z.png",res = 600, width=15, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Ly6k","Sycp3","Sycp1","Tex101","Emc7","Rad51ap2","Syce3","Sycp2","Meiob"))
dev.off()

png("./Spermatogenesis_DEGs/4.FeaturePlot_eP.png",res = 600, width=15, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Clgn", "Hormad1", "Cntd1", "Tdrd9", "Piwil2", "Scmh1", "Catsperb", "Kdm3a", "Msh4"))
dev.off()

png("./Spermatogenesis_DEGs/5.FeaturePlot_mP.png",res = 600, width=10, height=10, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Spata16","Spink2","Pgam2","Nme5"))
dev.off()

png("./Spermatogenesis_DEGs/6.FeaturePlot_lP.png",res = 600, width=10, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Ybx2","Pou5f2","Cdk1","Aurka","Spc25","Ccdc62"))
dev.off()

png("./Spermatogenesis_DEGs/7.FeaturePlot_RS.png",res = 600, width=15, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("Ssxb3","Ssxb1","Ssxb2","Zpbp","Acrv1", "Catsper3","Catsper4","Spata25","Izumo1"))
dev.off()

png("./Spermatogenesis_DEGs/8.FeaturePlot_ES.png",res = 600, width=15, height=15, units="in")
FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c("H1f9","H2ap","Prm1","Prm2","Prm3","Tnp1","Tnp2", "Gapdhs","Car2"))
dev.off()

new.cluster.ids <- c("L","P","pL","M","RS1","RS1","Z","RS1","P-like","D","Spg","pL","RS2")
names(new.cluster.ids) <- levels(Hsf5.merge.germ.cell)
Hsf5.merge.germ.cell <- RenameIdents(  Hsf5.merge.germ.cell, new.cluster.ids)

Hsf5.merge.germ.cell$celltype <- Hsf5.merge.germ.cell@active.ident
DimPlot(Hsf5.merge.germ.cell,reduction = "umap.harmony",label = TRUE,group.by = c("orig.ident", "celltype"))

# Manual annotation
# plot <- FeaturePlot(Hsf5.merge.germ.cell, reduction = "umap.harmony", features = c('Tubb5'))
# add_to_spg <- CellSelector(plot = plot)
# add_to_spg[! add_to_spg %in% names(as.list(Hsf5.merge.germ.cell@active.ident[Hsf5.merge.germ.cell@active.ident == "Spg"]))]
# ...

Hsf5.merge.germ.cell@active.ident <- factor(x=Hsf5.merge.germ.cell@active.ident, levels = c("Spg","pL","L","Z","eP","mP","P-like","lP","D","M","RS1","RS2"))
Hsf5.merge.germ.cell$celltype <- Hsf5.merge.germ.cell@active.ident

new.cluster.ids <- c("0", "1", "2", "3", "4","5", "6", "7", "8", "9", "10", "11")
names(new.cluster.ids) <- levels( Hsf5.merge.germ.cell)
Hsf5.merge.germ.cell <- RenameIdents(Hsf5.merge.germ.cell, new.cluster.ids)
Hsf5.merge.germ.cell$celltype_order <- Hsf5.merge.germ.cell@active.ident
Hsf5.merge.germ.cell@active.ident <- Hsf5.merge.germ.cell$celltype

png(file="./Germ_cells_identified.png",res = 600, width=15, height=15, units="in")
DimPlot(Hsf5.merge.germ.cell,
        reduction = "umap.harmony",
        group.by = c("celltype_order","celltype"),
        label = TRUE)
dev.off()
write_rds(Hsf5.merge.germ.cell, "../data/Germ_cells.rds")
DimPlot(Hsf5.merge.ident, reduction = "umap.harmony",label = TRUE)

## data visualization
Hsf5.merge.ident <- readRDS("../data/Whole_testis.rds")
Hsf5.merge.germ.cell <- readRDS("../data/Germ_cells.rds")

## Whole testis cell UMAP analysis 
p1 <- DimPlot(Hsf5.merge.ident, reduction = "umap.harmony",label = F) + xlim(-12, 12) + ylim(-13, 16) + labs(title = "") + xlab("UMAP_1") + ylab("UMAP_2") + NoLegend()
p2 <- DimPlot(subset(Hsf5.merge.ident, subset = orig.ident == "Hsf5.WT"), reduction = "umap.harmony", cols = c("Hsf5.KO" = "#f8766d", "Hsf5.WT" = "#00bfc4"), group.by = "orig.ident")
p3 <- DimPlot(subset(Hsf5.merge.ident, subset = orig.ident == "Hsf5.KO"), reduction = "umap.harmony", cols = c("Hsf5.KO" = "#f8766d", "Hsf5.WT" = "#00bfc4"), group.by = "orig.ident")

pdf(file = "../whole cell UMAP.pdf", width = 12, height = 4)
combined_plot <- grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

## Whole testis marker gene heatmap
markers <- FindAllMarkers(JoinLayers(Hsf5.merge.ident), only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = 1)
write.csv(markers,"../figS1/data/markers_all.csv")

Marker <- markers
###Calculate Average Expression Value
FS <- names(table(Marker$gene)[table(Marker$gene)==1])
Marker.f <- subset(Marker,gene %in% FS)
FS_list <- NULL;for(i in c("SPG","SCytes","STids",'Sertoli cells','Leydig cells', 'Macrophages','Endothelial cells','Interstitial progenitors', 'Peritubular myoid cells',  'Perivascular cells')){
    genes <- rownames(subset(Marker.f, cluster %in% i))
    if (length(genes) > 200) {
        genes <- head(genes, 200)
    }
    FS_list <- c(FS_list, genes)
}

Exp <- AggregateExpression(JoinLayers(Hsf5.merge.ident), assays = "RNA")
Exp <- Exp$RNA
Exp <- Exp[FS_list,]
Exp <- Exp[,c("SPG","SCytes","STids",'Sertoli cells','Leydig cells', 'Macrophages','Endothelial cells','Interstitial progenitors', 'Peritubular myoid cells',  'Perivascular cells')]

pdf(file="../Whole testis marker gene.pdf", width = 4, height = 6)
pheatmap(Exp,cluster_rows = FALSE,cluster_cols = FALSE,scale="row",show_rownames=FALSE, angle_col = 45)
dev.off()


## Sperm cell marker gene violin plot
library(MySeuratWrappers)
markers <- FindAllMarkers(JoinLayers(Hsf5.merge.germ.cell), only.pos = TRUE, min.diff.pct = 0.2, logfc.threshold = 1)
write.csv(markers,"../data/germ_cells_markers_all.csv")

Vln_markers <- c( 'Ccnd2',"Dmrt1", 'Prdm9', 'Rad51ap2',"Mnd1","Mcmdc2","H1f6", "Pomc", 'Pou5f2',"Rhcg", 'Tex36', 'Sun5', 'Cstl1')

pdf("../Sperm cell marker gene violin plot.pdf", width=8, height=2.5)
VlnPlot(JoinLayers(subset(Hsf5.merge.germ.cell, subset = celltype %in% c("Spg","pL","L","Z","eP","mP","lP","D","M","RS1","RS2"))),
        # JoinLayers(subset(markerdata, subset = orig.ident == "Hsf5.WT")),
        features = Vln_markers,
        stacked=T,
        pt.size= 0,
        direction = "horizontal", #水平作图
        x.lab = '', y.lab = '')
dev.off()


## Whole testis marker gene FeaturePlot
aspect_ratio = 0.8

## Sertoli
p1 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Sox9"), order = T, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Leydig cells
p2 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Star"), order = T, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Interstitial progenitors
p3 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Tcf21"), order = T, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Macrophages 10.3389/fimmu.2018.02246
p4 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Adgre1"), order = T, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Endothelial cells
p5 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Egfl7"), order = T, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Peritubular myoid cells PTMs
p6 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Acta2"), order = F, aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Perivascular cells 
p7 <- FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Angpt1"), order = T , aspect_ratio = aspect_ratio, cols = c("lightgrey", "red")) + xlab("UMAP_1") + ylab("UMAP_2") + theme(legend.position = "none")

## Legend 
legend <- cowplot::get_legend(FeaturePlot_scCustom(Hsf5.merge.ident, reduction = "umap.harmony", features = c("Angpt1")) + Move_Legend("right"))

pdf(file = "../Whole testis marker gene FeaturePlot.pdf", width = 12, height = 5.5)
combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, legend, ncol = 4)
dev.off()


## Remapping analysis of reference cell annotation sets
file_paths <- list.files("../data/chen/GSE107644_RAW", pattern = "\\.txt$", full.names = TRUE)
data <- read.table(file_paths[1], header = TRUE, row.names = 1, sep = "\t")
prefix <- strsplit(gsub(".txt", "", basename(file_paths[1])),"_" )[[1]][2] # Extract prefix from file name
colnames(data) <- gsub(paste0("^", prefix, "_"), "", colnames(data))
chen_reference <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, assay = "RNA")

# Loop through each count matrix file and add it to a Seurat object
for (file_path in file_paths[-1]) {
    data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
    prefix <- strsplit(gsub(".txt", "", basename(file_path)),"_" )[[1]][2] # Extract prefix from file name
    colnames(data) <- gsub(paste0("^", prefix, "_"), "", colnames(data))
    chen_reference <- merge(chen_reference, CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, assay = "RNA"))
}

chen_reference <- JoinLayers(chen_reference)

chen_reference <- NormalizeData(chen_reference, normalization.method = "LogNormalize", scale.factor = 10000)
chen_reference <- FindVariableFeatures(chen_reference, selection.method = "vst", nfeatures = 2000)
chen_reference <- ScaleData(chen_reference)
chen_reference <- RunPCA(chen_reference, features = VariableFeatures(object = chen_reference))
chen_reference <- FindNeighbors(chen_reference, dims = 1:16, reduction = "pca")
chen_reference <- FindClusters(chen_reference, resolution = 0.5, cluster.name = "unintegrated_clusters")
chen_reference <- RunUMAP(chen_reference, dims = 1:16, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(chen_reference, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

anchors <- FindTransferAnchors(
    reference = Hsf5.merge.germ.cell,
    query = chen_reference,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    dims = 1:16
)

chen_reference <- MapQuery(
    anchorset = anchors,
    query = chen_reference,
    reference = Hsf5.merge.germ.cell,
    refdata = list(
        celltype.1 = "celltype"
    ),
    reference.reduction = "pca", 
    reduction.model = "umap.harmony",
    verbose = TRUE
)

p1 <- DimPlot(Hsf5.merge.germ.cell,
            reduction = "umap.harmony",
            group.by = c("celltype"),
            cols =  c("Spg"="#52a95f","pL"="#56c4f3","L"="#7BAFDE","Z"="#9db7a5","eP"="#9b77a2","mP"="#bbbbff","P-like"="#ff6e5f","lP"="#fbb1a2","D"="#e69dc5","M"="#a6d608","RS1"="#F2A874","RS2"="#456d88"),
            label = F,
            label.box = F,
            label.size = 3,
            repel = TRUE) + NoLegend() + xlab("UMAP_1") + ylab("UMAP_2") + xlim(-11, 10) + ylim(-11, 14) + labs(title = "") #  + ggtitle("Hsf5's predict cell")
p2 <- Seurat::DimPlot(chen_reference,
            reduction = "ref.umap",
            group.by = c("orig.ident"),
            pt.size = 2,
            label = T,
            label.box = T,
            label.size = 4.25,
            label.color = "white",
            repel = TRUE) + xlab("UMAP_1") + ylab("UMAP_2") + NoLegend() + xlim(-11, 10) + ylim(-11, 14) + labs(title = "")

pdf(file="./Reference gene set remapping.pdf", width=12, height=6)
grid.arrange(p1, p2, ncol = 2)
dev.off()

## Cell type count
library(ggbreak)
# Count the number of cells from different samples in each cluster
cluster_table <- table(Hsf5.merge.germ.cell@active.ident, Hsf5.merge.germ.cell$orig.ident)

cluster_df <- as.data.frame.matrix(cluster_table)
colnames(cluster_df) <- c("Hsf5", "Wildtype")
cluster_df$celltype = row.names(cluster_df)

cell_types <- c("Spg","pL","L","Z","eP","mP","lP","D","M","RS1","RS2","P-like")  
fill_colors <- c("Spg"="#52a95f","pL"="#56c4f3","L"="#7BAFDE","Z"="#9db7a5","eP"="#9b77a2","mP"="#bbbbff","P-like"="#ff6e5f","lP"="#fbb1a2","D"="#e69dc5","M"="#a6d608","RS1"="#F2A874","RS2"="#456d88")

df = data.frame(sample = c(rep("Wildtype",12),rep("Hsf5-/-",12)))
df$celltype <- factor(c(cluster_df$celltype, cluster_df$celltype),levels = rev(cell_types))
df$counts <- c(cluster_df$Wildtype, -cluster_df$Hsf5)

pdf(file="./Cell type count.pdf", width=4, height=4)
ggplot(df, aes(x=counts, y=celltype, fill=celltype)) +
    scale_fill_manual(values = fill_colors) +
    geom_bar(stat = "identity") +theme_classic() +
    theme(title = element_blank() , axis.text.y = element_text(color = "black",size = 12) , axis.text.x = element_text(color = "black",size = 12)) + 
    NoLegend() +
    geom_text(aes(label=counts), hjust=-0.2, vjust=0.5, size=4, color="black", position=position_dodge(width=0.9)) +
    coord_cartesian(xlim = c(-700, 700)) + 
    geom_vline(xintercept = 0, linetype = "solid")
dev.off()

## Export data for RNA velocity analysis
library(monocle3)
library(SeuratWrappers)

df<-data.frame(Cells=Cells(Hsf5.merge.germ.cell), sample=Hsf5.merge.germ.cell$sample) 
df$id<-sapply(df$Cells,function(x)paste(unlist(strsplit(x, "-"))[1],"x",sep = ""))
df$Cells<- gsub("_", ":", df$id)
write.csv(df$Cells, file = "./velocyto/cell_ID_obs.csv", row.names = FALSE) 

cell_embeddings<-Embeddings(Hsf5.merge.germ.cell, reduction = "umap.harmony")
rownames(cell_embeddings)<-df$Cells 
write.csv(cell_embeddings, file = "./velocyto/cell_embeddings.csv")

clusters_obs<-Hsf5.merge.germ.cell$celltype 
names(clusters_obs)<-df$Cells 
write.csv(clusters_obs, file = "./velocyto/cell_type_obs.csv")

clusters_obs<-Hsf5.merge.germ.cell$celltype_order 
names(clusters_obs)<-df$Cells 
write.csv(clusters_obs, file = "./velocyto/cell_clusters_obs.csv")


## Analysis of Pachytene stage
Hsf5.merge.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("eP","mP","lP","P-like"))
Hsf5.merge.pachytene <- FindVariableFeatures(Hsf5.merge.pachytene, selection.method = "vst", nfeatures = 1000)
Hsf5.merge.pachytene <- RunPCA(Hsf5.merge.pachytene, features = VariableFeatures(object = Hsf5.merge.pachytene)[! VariableFeatures(object = Hsf5.merge.pachytene) %in% c("Ftl1")])
Hsf5.merge.pachytene <- IntegrateLayers(object = Hsf5.merge.pachytene, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",verbose = FALSE)
Hsf5.merge.pachytene <- FindNeighbors(Hsf5.merge.pachytene, dims = 1:4, reduction = "harmony")
Hsf5.merge.pachytene <- FindClusters(Hsf5.merge.pachytene, resolution = 0.5, cluster.name = "harmony_clusters")
Hsf5.merge.pachytene <- RunUMAP(Hsf5.merge.pachytene, dims = 1:4, reduction = "harmony", reduction.name = "umap.harmony")
DimPlot(Hsf5.merge.pachytene, reduction = "umap.harmony", group.by = c("orig.ident", "celltype"))

Hsf5.merge.pachytene@active.ident <- Hsf5.merge.pachytene$celltype

DimPlot(Hsf5.merge.pachytene, reduction = "umap.harmony", group.by = c("orig.ident", "celltype"))
write_rds(Hsf5.merge.pachytene, "./data/Pachytene_cells.rds")

p1 <- Seurat::DimPlot(Hsf5.merge.pachytene, reduction = "umap.harmony", group.by = "celltype", cols =  c("eP"="#9b77a2","mP"="#bbbbff","P-like"="#ff6e5f","lP"="#fbb1a2"),label = T,label.box = T,label.color = "white") + xlim(-14, 10) + ylim(-6, 6) + labs(title = "") + xlab("UMAP_1") + ylab("UMAP_2") + NoLegend()
p2 <- Seurat::DimPlot(Hsf5.merge.pachytene, reduction = "umap.harmony", cols = c("Hsf5.KO" = "#f8766d", "Hsf5.WT" = "#00bfc4"),label = T,label.box = T,label.color = "white", group.by = "orig.ident") + NoLegend() + xlim(-14, 10) + ylim(-6, 6) + labs(title = "") + xlab("UMAP_1") + ylab("UMAP_2")

pdf(file = "../Pachytene_UMAP.pdf", width = 7, height = 3.5)
combined_plot <- grid.arrange(p1, p2, ncol = 2)
dev.off()

Hsf5.merge.pachytene <- readRDS("./data/Pachytene_cells.rds")

Hsf5.wt.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("eP", "mP", "lP") & orig.ident == "Hsf5.WT")
Hsf5.ko.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("P-like") & orig.ident == "Hsf5.KO")

Hsf5.merge.pachytene <- merge(Hsf5.wt.pachytene, y = Hsf5.ko.pachytene)

## Differential Gene Volcano Plot
P_like_DEG <- FindMarkers(JoinLayers(Hsf5.merge.pachytene), ident.1 = "P-like", ident.2 = "mP",test.use = "MAST")
P_like_DEG <- P_like_DEG %>%
    mutate(Difference = pct.1 - pct.2) %>% 
        rownames_to_column("gene")

log2FC = 1
padj = 0.01

P_like_DEG$threshold="ns";
P_like_DEG[which(P_like_DEG$avg_log2FC  > log2FC & P_like_DEG$p_val_adj <padj),]$threshold="up";
P_like_DEG[which(P_like_DEG$avg_log2FC  < (-log2FC) & P_like_DEG$p_val_adj < padj),]$threshold="down";
P_like_DEG$threshold=factor(P_like_DEG$threshold, levels=c('down','ns','up'))

P_like_DEG <- P_like_DEG[order(P_like_DEG$p_val_adj, decreasing = TRUE), ]

p <- ggplot(P_like_DEG, aes(x=Difference, y=avg_log2FC,color = threshold)) + 
    geom_point(size=1, alpha = 0.3) +
    geom_point(color = "red", size = 1.5, data = subset(P_like_DEG, gene %in% c("Hsf2","Hspa12b", "Dmc1","Ddx11","Ubr2","Msh4","Meiob","Hat1","Cenpp","Tex15","Spata2","Scmh1","Ctsl","Clock","Smad4","Piwil2")))+
    geom_point(color = "blue",size = 1.5 ,data = subset(P_like_DEG, gene %in% c("Hspb9","Hspa1l","Hspb6","Hsp90aa1","Hspb11", "Hspa2","Upf3a","H2ax","Ccna1","Mns1","Krt88")))+
    geom_vline(xintercept = 0.0,linetype=2)+
    geom_hline(yintercept = 0,linetype=2)+
    scale_color_manual(values=c("#5fb9ff","grey","#de7a60"))+
    geom_label_repel(data=subset(P_like_DEG, gene %in% c("Hsf2","Hspa12b","Dmc1","Ddx11","Ubr2","Msh4","Meiob","Hat1","Cenpp","Hsf2","Tex15","Spata2","Scmh1","Ctsl","Clock","Smad4","Piwil2")), 
                    aes(label=gene),  
                    color="#b20000",
                    label.padding = 0.1,
                    box.padding = 0.35,
                    max.overlaps = Inf,
                    segment.size = 0.4,
                    min.segment.length = 0,
                    size=4)+
    geom_label_repel(data=subset(P_like_DEG, gene %in% c("Hspb9","Hspa1l","Hspb6","Hsp90aa1","Hspb11", "Hspa2","Upf3a","H2ax","Ccna1","Mns1","Krt88")), 
                    aes(label=gene), 
                    color="#0000b2", 
                    label.padding = 0.1,
                    box.padding = 0.35,
                    max.overlaps = Inf,
                    segment.size = 0.4,
                    min.segment.length = 0,
                    size=4)+theme_classic()

pdf(file = "./DEG_of_KO-P-like_vs_WT_mP.pdf", width = 4.5, height = 5.5)
p
dev.off()


## Venn diagram of differential genes
Hsf5.wt.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("eP", "mP", "lP") & orig.ident == "Hsf5.WT")
Hsf5.ko.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("eP","mP", "P-like") & orig.ident == "Hsf5.KO")

P_like_vs_eP_DEG <- FindMarkers(Hsf5.ko.pachytene, ident.1 = "P-like", ident.2 = "eP",test.use = "MAST")
mP_vs_eP_DEG <- FindMarkers(Hsf5.wt.pachytene, ident.1 = "mP", ident.2 = "eP",test.use = "MAST")
log2FC = 1
write.csv(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC)),"./Hsf5-KO_P_like_vs_eP_Up-regulated_DEGs.csv")
write.csv(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC)),"./WT_mP_vs_eP_Up-regulated_DEGs.csv.csv")
write.csv(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC)),"./Hsf5-KO_P_like_vs_eP_Down-regulated_DEGs.csv")
write.csv(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC)),"./WT_mP_vs_eP_Down-regulated_DEGs.csv.csv")

Hsf5.wt.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("eP", "mP", "lP") & orig.ident == "Hsf5.WT")
Hsf5.ko.pachytene <- subset(Hsf5.merge.germ.cell, subset = celltype %in% c("P-like") & orig.ident == "Hsf5.KO")

Hsf5.merge.pachytene <- merge(Hsf5.wt.pachytene, y = Hsf5.ko.pachytene)

P_like_DEG <- FindMarkers(JoinLayers(Hsf5.merge.pachytene), ident.1 = "P-like", ident.2 = "mP",test.use = "MAST")
write.csv(subset(P_like_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC)),file = "./Hsf5-KO_P_like_vs_WT_mP_Up-regulated_DEGs.csv")
write.csv(subset(P_like_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC)),file = "./Hsf5-KO_P_like_vs_WT_mP_Down-regulated_DEGs.csv")

log2FC = 1

max_length <- max(length(KO_Plike_up <- rownames(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC)))),
            length(KO_Plike_down <- rownames(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC)))),
            length(WT_mP_up <- rownames(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC)))),
            length(WT_mP_down <- rownames(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC)))))

DEG <- list("Hsf5-/- P-like vs eP up-regulated"=rownames(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC))),
            "Hsf5-/- P-like vs eP down-regulated"=rownames(subset(P_like_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC))),
            "Hsf5+/+ mP vs eP up-regulated"=rownames(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC > log2FC))),
            "Hsf5+/+ mP vs eP down-regulated"=rownames(subset(mP_vs_eP_DEG, p_val_adj<0.01 & (avg_log2FC < -log2FC))))

library(ggvenn)
ggvenn(
    data = DEG,
    columns = c("Hsf5-/- P-like vs eP up-regulated","Hsf5+/+ mP vs eP up-regulated"),
    show_elements = F,
    label_sep = "\n",
    show_percentage = F,
    fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"),
    fill_alpha = 0.5,
    stroke_color = "white",
    stroke_alpha = 0.5,
    stroke_size = 0.5
    stroke_linetype = "solid",
    set_name_color = "black",
    set_name_size = 6,
    text_color = "black",
    text_size = 4
)


ggvenn(
    data = DEG,
    columns = c("Hsf5-/- P-like vs eP down-regulated","Hsf5+/+ mP vs eP down-regulated"),
    show_elements = F,
    label_sep = "\n",
    show_percentage = F,
    fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"),
    fill_alpha = 0.5,
    stroke_color = "white",
    stroke_alpha = 0.5,
    stroke_size = 0.5
    stroke_linetype = "solid",
    set_name_color = "black",
    set_name_size = 6,
    text_color = "black",
    text_size = 4
)

library(eulerr)
plot(euler(c("P-like"=509,"mP"=311,"P-like&mP"=69),
    shape = "ellipse"),
    fills = list(fill=c("#E41A1C", "#1E90FF", "#646993", "#80FF00"),
    alpha=0.5),
    quantities = c(509,311,69,col="black",font=1,cex=4),
    labels = list(col="white",font=1,cex=4))

pdf(file = "./eulerr_up.pdf", width = 3, height = 2)
plot(euler(c("P-like"=663,"mP"=257,"P-like&mP"=159),
    shape = "ellipse"),
    fills = list(fill=c("#E41A1C", "#1E90FF", "#646993", "#80FF00"),
    alpha=0.5),
    quantities = c(663,257,159,col="black",font=1,cex=4))
dev.off()

pdf(file = "./eulerr_down.pdf", width = 3, height = 2)
plot(euler(c("P-like"=509,"mP"=311,"P-like&mP"=69),
    shape = "ellipse"),
    fills = list(fill=c("#E41A1C", "#1E90FF", "#646993", "#80FF00"),
    alpha=0.5),
    quantities = c(509,311,69,col="black",font=1,cex=4))
dev.off()

## Export data for RNA velocity analysis
library(monocle3)
library(SeuratWrappers)

Hsf5.merge.pachytene <- subset(Hsf5.merge.pachytene, subset = celltype %in% c("eP","mP","lP","P-like"))

df<-data.frame(Cells=Cells(Hsf5.merge.pachytene), sample=Hsf5.merge.pachytene$sample) 
df$id<-sapply(df$Cells,function(x)paste(unlist(strsplit(x, "-"))[1],"x",sep = ""))
df$Cells<- gsub("_", ":", df$id)
write.csv(df$Cells, file = "./velocyto/pachytene_cell_ID_obs.csv", row.names = FALSE) 

cell_embeddings<-Embeddings(Hsf5.merge.pachytene, reduction = "umap.harmony")
rownames(cell_embeddings)<-df$Cells 
write.csv(cell_embeddings, file = "./velocyto/pachytene_cell_embeddings.csv")

clusters_obs<-Hsf5.merge.pachytene$celltype 
names(clusters_obs)<-df$Cells 
write.csv(clusters_obs, file = "./velocyto/pachytene_cell_type_obs.csv")

clusters_obs<-Hsf5.merge.pachytene$celltype_order 
names(clusters_obs)<-df$Cells 
write.csv(clusters_obs, file = "./velocyto/pachytene_cell_clusters_obs.csv")