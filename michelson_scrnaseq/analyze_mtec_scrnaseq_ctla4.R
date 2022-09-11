#####
#load data
library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)
library(ggthemes)
library(BuenColors)
library(ggrepel)

setwd( "C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_100821" )
seurat_table <- readRDS( "Rdata/meclo2_seurat_table_2021-11-23.rds" )

date="20220328"

cmap <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]$`Tableau 20`$value

cmap2 = jdb_color_maps
cmap2 <- as.vector(cmap2)
cmap2 <- c(cmap2,"darkseagreen4")

names(cmap2) <- c("Aire-expressing","Tuft","Microfold","Enterocyte/hepatocyte","Neuroendocrine","Immature","Ionocyte","Aire-deficient","Ptf1a+ ductal",
                  "Keratinocyte","TA","Basal","Ciliated","Lung","perinatal cTEC","Muscle","adult cTEC","Tuft2")

cmap2["adult cTEC"] <- "darkorchid4"
cmap2["Lung"] <- "#993299"
cmap2["Neuroendocrine"] <- "#79adba"
cmap2["Immature"] <- "#f9838b"

cmap <- cmap2

seurat_table <- seurat_table[,seurat_table$geno.ident %in% c("wt","aire")]

#####
#viz data

feature <- "Irf8"

a <- DimPlot(seurat_table, reduction="umap", label=T, cols=cmap, repel=T ) + theme(legend.position="none")
x <- FeaturePlot(seurat_table, features = c(feature), pt.size=1, order=T ) + theme(legend.position="none")
y <- VlnPlot( seurat_table, features=c(feature), cols=cmap ) + theme(legend.position="none")

pdf( paste0("figures/umap_",feature, "_",date,".pdf"), height=4, width=5 )
plot_grid(a,x,y, ncol=3)
# x
dev.off()

pdf( paste0("figures/featureScatter_Ctla4_Irf8_",date,".pdf"), height=4, width=6 )
FeatureScatter(seurat_table, feature1="Ctla4",feature2="Irf8", cols=cmap)
dev.off()

#####
#analyze cluster composition
cluster.markers <- FindMarkers(seurat_table, ident.1="Aire-deficient", min.pct=0.1, logfc.threshold=1, only.pos=F, test.use="bimod")
head(cluster.markers, n = 20)

#####
#analyze cluster composition

idx1 <- seurat_table@assays$RNA@counts["Ctla4",]>0
idx2 <- sample(seurat_table@assays$RNA@counts["Ctla4",]==0, size=sum(seurat_table@assays$RNA@counts["Ctla4",]>0))

set.seed(12345)

cluster.markers.ctla4 <- FindMarkers(seurat_table@assays$RNA, 
                                    cells.1=colnames(seurat_table)[idx1], 
                                    cells.2=names(idx2), 
                                    only.pos=F, logfc.threshold=0.0, min.pct=0, test.use="DESeq2")

cluster.markers <- cluster.markers.ctla4
label_idx <- ifelse( abs(cluster.markers$avg_log2FC) > 1 & -log10(cluster.markers$p_val_adj) > 1.3, rownames(cluster.markers), "" )
color_idx <- ifelse(rownames(cluster.markers)!="" & cluster.markers$avg_log2FC > 1 & -log10(cluster.markers$p_val_adj) > 1.3,"1","0")
color_idx <- ifelse(rownames(cluster.markers)!="" & cluster.markers$avg_log2FC < -1 & -log10(cluster.markers$p_val_adj) > 1.3,"2",color_idx)

pdf("figures/ctla4_vs_subsampledAll_DESeq2_volcano_20220328.pdf", height=4, width=6.5)
ggplot( cluster.markers, aes(x=avg_log2FC, y=-log10(p_val_adj) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel(label=label_idx ) +
  scale_color_manual( values=c("gray","firebrick3","navyblue") ) +
  geom_hline( yintercept=1.3, lty="dashed" ) +
  geom_vline( xintercept=1, lty="dashed" ) +
  geom_vline( xintercept=-1, lty="dashed" ) +
  theme_few() +
  scale_x_continuous(limits=c(-3,3), oob=squish) +
  scale_y_continuous(limits=c(0,10), oob=squish) +
  theme(legend.position="none") +
  xlab( paste0("log2 FC") ) +
  ylab( paste0("-log10 adjusted p-value") ) +
  ggtitle("")
dev.off()

write_delim( as.data.frame(rownames(cluster.markers)[cluster.markers$avg_log2FC > 1 & -log10(cluster.markers$p_val_adj) > 1.3]), "text_outputs/ctla4_UP_genes_20220328.txt", col_names=F, delim="\t" )
write_delim( as.data.frame(rownames(cluster.markers)[cluster.markers$avg_log2FC < -1 & -log10(cluster.markers$p_val_adj) > 1.3]), "text_outputs/ctla4_DN_genes_20220328.txt", col_names=F, delim="\t" )

gprofiler <- read.delim( "text_outputs/20220328_ctla4/gProfiler_mmusculus_3-28-2022_11-25-49 AM__intersections.csv", header=T, sep="," )

head(gprofiler)
gprofiler <- gprofiler[gprofiler[,1]=="GO:BP",]
gprofiler$term_name <- factor(gprofiler$term_name, levels=rev(c(gprofiler$term_name)))

pdf("figures/ctla4_vs_subsampledAll_DESeq2_pathway_full_20220328.pdf", height=12, width=9)
ggplot(gprofiler, aes(x=negative_log10_of_adjusted_p_value, y=term_name)) +
  geom_col(color="black", fill="gray") +
  theme_few()
dev.off()
pdf("figures/ctla4_vs_subsampledAll_DESeq2_pathway_select_20220328.pdf", height=2.75, width=9)
ggplot(gprofiler[c(1,14,16,17,21),], aes(x=negative_log10_of_adjusted_p_value, y=term_name)) +
  geom_col(color="black", fill="gray") +
  theme_few()
dev.off()
