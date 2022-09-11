# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Aug 15 12:52:01 EDT 2019

setwd( "C:/Users/dmich/Dropbox (Personal)/CBDM Lab/Data/MEC_Ctla4/Microarray_and_NOD" )

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE8564", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "22233300001111444555"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G5-G0, G1-G0, G2-G1, G3-G2, G4-G3, G5-G4, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
# write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE8564", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1261", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "22233300001111444555"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("BALB+WT","BALB+KO","B6+WT","B6+KO","NOD+WT","NOD+KO")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#dcdaa5","#f2cb98","#dff4e4","#f4dff4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE8564", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")

#####################################################
#plot Ctla4 data

library(tidyverse)
library(ggplot2)

ctla4_exp <- as.numeric( exprs(gset)["1419334_at",1:20] )
sample_id <- c("B6 WT1","B6 WT2","B6 WT3","B6 KO1","B6 KO2","B6 KO3",
               "BALB WT1","BALB WT2","BALB WT3","BALB WT4","BALB KO1","BALB KO2","BALB KO3","BALB KO4",
               "NOD WT1", "NOD WT2","NOD WT3","NOD KO1","NOD KO2","NOD KO3")
group_id <- factor(c("B6\nWT","B6\nWT","B6\nWT","B6\nKO","B6\nKO","B6\nKO",
               "BALB\nWT","BALB\nWT","BALB\nWT","BALB\nWT","BALB\nKO","BALB\nKO","BALB\nKO","BALB\nKO",
               "NOD\nWT", "NOD\nWT","NOD\nWT","NOD\nKO","NOD\nKO","NOD\nKO"), levels=c("B6\nWT","B6\nKO",
                                                                                       "BALB\nWT","BALB\nKO",
                                                                                       "NOD\nWT","NOD\nKO"), ordered=T )
strain_id <- factor(c("C57BL/6","C57BL/6","C57BL/6","C57BL/6","C57BL/6","C57BL/6",
              "BALB/c","BALB/c","BALB/c","BALB/c","BALB/c","BALB/c","BALB/c","BALB/c",
              "NOD", "NOD","NOD","NOD","NOD","NOD"), levels=c("C57BL/6","BALB/c","NOD"), ordered=T)
ctla4_dataframe <- data.frame( sample_id = sample_id, 
                               group_id = as.factor( group_id ),
                               ctla4_exp =ctla4_exp,
                               strain_id= as.factor(strain_id))


ctla4_plot <- ggplot( ctla4_dataframe, aes( x=group_id, y=ctla4_exp, fill=group_id) ) +
  # geom_bar( stat="identity" ) +
  geom_boxplot() +
  geom_jitter() +
  xlab( "" ) +
  ylab( "Ctla4 expression (AU)") +
  coord_trans( y="log2") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0), legend.position = "none") +
  scale_fill_manual( values=c("steelblue4","steelblue1","grey28","grey55","plum4","plum2") )

pdf( "20220827_microarray_ctla4_expression_perStrain.pdf", height=4, width=6 )
ctla4_plot + facet_wrap(.~strain_id, scales="free")
dev.off()

#analyze expression fc

fc_b6 <- ctla4_dataframe[ctla4_dataframe[,"group_id"]=="B6\nKO","ctla4_exp"] /
  ctla4_dataframe[ctla4_dataframe[,"group_id"]=="B6\nWT","ctla4_exp"]

fc_balb <- ctla4_dataframe[ctla4_dataframe[,"group_id"]=="BALB\nKO","ctla4_exp"] /
  ctla4_dataframe[ctla4_dataframe[,"group_id"]=="BALB\nWT","ctla4_exp"]

fc_nod <- ctla4_dataframe[ctla4_dataframe[,"group_id"]=="NOD\nKO","ctla4_exp"] /
  ctla4_dataframe[ctla4_dataframe[,"group_id"]=="NOD\nWT","ctla4_exp"]

mean(fc_b6)
mean(fc_balb)
mean(fc_nod)

t.test( log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="B6\nKO","ctla4_exp"]),
        log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="B6\nWT","ctla4_exp"]),
        paired = F, alternative = "two.sided" )

t.test( log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="BALB\nKO","ctla4_exp"]),
        log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="BALB\nWT","ctla4_exp"]),
        paired = F, alternative = "two.sided" )

t.test( log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="NOD\nKO","ctla4_exp"]),
        log2(ctla4_dataframe[ctla4_dataframe[,"group_id"]=="NOD\nWT","ctla4_exp"]),
        paired = F, alternative = "two.sided" )
