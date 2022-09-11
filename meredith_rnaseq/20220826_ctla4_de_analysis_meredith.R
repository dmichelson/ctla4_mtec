library(tidyverse)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(pcaMethods)
library(ggrepel)
library(scales)
library(ggthemes)

setwd("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Scratch/20220826_ctla4/stag2_rnaseq")
date = "20220826"

count_table <- read.delim( "GSE180935_Gene_count_table.tsv", sep="\t", header=T, row.names=1 )
count_table <- count_table[,c(1,2,3,8,9,10)]

colnames(count_table) <- c("AireWT1","AireWT2","AireWT3","AireKO1","AireKO2","AireKO3")
sample_names <- colnames(count_table)

#####
##de analysis

# run edgeR pipeline
i = "AireWT"
j = "AireKO"

idx <- substr(colnames(count_table),1,6)== i
ctrl_idx <- substr(colnames(count_table),1,6)== j

group <- factor(c(rep(1,times=sum(idx)),rep(2,times=sum(ctrl_idx))))
y <- DGEList( counts=cbind( count_table[,idx], count_table[,ctrl_idx]), 
              group=group, 
              genes=rownames(count_table) )

keep <- filterByExpr(y, min.count=10)
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

design <- model.matrix(~group)
y <- estimateDisp(y, design)

fit <- glmQLFit(y,design)

qlf <- glmQLFTest(fit,coef=2)

qlf$table <- cbind( qlf$table, p.adjust( qlf$table$PValue, method="BH" ) )
colnames(qlf$table)[5] = "q_val"

write_delim( cbind( rownames(qlf$table),qlf$table ), paste0("de_analysis/",i,"_vs_", j,"_DE_table_", date ,".txt"), col_names=T, delim="\t" )


#plot volcano
label_idx <- ifelse( rownames(qlf$table) %in% c("Ctla4"), rownames(qlf$table), "" )
color_idx <- ifelse( rownames(qlf$table) %in% rownames(topTags(qlf, n=25)), "1", "0" )


pdf(paste0("figures/",i,"_vs_", j,"_volcano_unlab_", date ,".pdf"), height=6, width=8)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( ) +
  # geom_point( data=qlf$table[features,], aes(x=qlf$table[features,"logFC"], y=-log10(qlf$table[features,"PValue"])), color="red" ) +
  # geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  # scale_x_continuous(limits=c(-6,6), oob=squish) +
  # scale_y_continuous(limits=c(0,8), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

pdf(paste0("figures/",i,"_vs_", j,"_volcano_lab_", date ,".pdf"), height=6, width=8)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$PValue) ) ) +
  geom_point( aes(color=color_idx) ) +
  # geom_point( data=qlf$table[features,], aes(x=qlf$table[features,"logFC"], y=-log10(qlf$table[features,"PValue"])), color="red" ) +
  geom_text_repel( label=label_idx, color="red", segment.alpha=0.75, segment.color="gray", segment.size=0.25 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 P-value") +
  scale_color_manual( values=c("black","red") ) +
  theme_bw() +
  # scale_x_continuous(limits=c(-6,6), oob=squish) +
  # scale_y_continuous(limits=c(0,8), oob=squish) +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
dev.off()

fc_cutoff = 2
qval_cutoff = 0.05

de_genes <- abs( qlf$table$logFC ) > fc_cutoff & qlf$table$q_val < qval_cutoff
down_sig <- qlf$table$logFC > fc_cutoff & qlf$table$q_val < qval_cutoff
up_sig <- qlf$table$logFC < -fc_cutoff & qlf$table$q_val < qval_cutoff

label_idx <- ifelse( rownames(qlf$table) %in% c("Ctla4"), rownames(qlf$table), "" )
color_idx <- vector( length=nrow(qlf$table) )
for (k in 1:nrow(qlf$table)) {
  if (rownames(qlf$table)[k] %in% rownames(qlf$table)[up_sig]){
    color_idx[k] <- "up"
  } else if (rownames(qlf$table)[k] %in% rownames(qlf$table)[down_sig]) {
    color_idx[k] <- "down"
  } else {
    color_idx[k] <- "neutral"
  }
}

pdf(paste0("figures/",i,"_vs_", j,"_volcano_up_down_sig_logfc",fc_cutoff,"_qval",qval_cutoff,"_", date ,".pdf"), height=4, width=6)
ggplot( qlf$table, aes( x=-qlf$table$logFC, y=-log10(qlf$table$q_val) ) ) +
  geom_point( aes(color=color_idx) ) +
  geom_text_repel( label=label_idx, segment.alpha=1, segment.color="black", segment.size=0.25, nudge_x = -1, nudge_y = 0.5 ) +
  xlab(paste0("log2 fold change, ", i, " vs ", j)) +
  ylab("-log10 adjusted P-value") +
  scale_color_manual( values=c("blue","black","red") ) +
  theme_bw() +
  theme( legend.position="none" ) +
  ggtitle( paste0(i, " vs ", j) )
  # scale_x_continuous(limits=c(-6,6), oob=squish) +
  # scale_y_continuous(limits=c(0,8), oob=squish) +
  # annotate( geom="text", label=sum(up_sig), x=6, y=7, color="blue" ) +
  # annotate( geom="text", label=sum(down_sig), x=-6, y=7, color="red" )
dev.off()
