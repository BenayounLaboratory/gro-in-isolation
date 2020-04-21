##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-04-22

#######        Week 5 code         #######

### Example of analysis of a microarray experiment  (Part 1) ###

###################################################################################################################
# Dataset for Example: Calorie restriction in different organs
# Publicly available from GEO Datasets repository
#### http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11291
#### Series GSE11291 (Public on Jun 10, 2008)
#### Effect of age, calorie restriction and resveratrol on gene expression in mouse heart, brain, and skeletal muscle
#### Platform: GPL1261 [Mouse430_2] Affymetrix Mouse Genome 430 2.0 Array
###################################################################################################################

###################################################################################################################
### Set the working directory
setwd('/Users/berenice/Dropbox/USC_Work_Folder/Teaching/2019-2020/GRO-in-isolation/Github/Week5')
options(stringsAsFactors = F)

###################################################################################################################
### 1. RMA (Robust Microarray Analysis) normalization of affy data
###################################################################################################################
library('affy')

# Read files in from directory that have relevant extension
files <- paste0("GSE11291_RAW/",list.files(path = "GSE11291_RAW/", pattern = "\\.CEL$"))

# RMA is a background correction and cross-sample normalization method for microarray
# just.rma will perform RMA normalization and load probe annotation from web using CEL file metadata
CR.Data <- just.rma(files)

# Write normalized data to external file, write an ExpressionSet to file
my.output <- paste0(Sys.Date(),"_GSE11291_CR_Data.RMA.Data.txt")
write.exprs(CR.Data, file = my.output, sep = "\t", eol="\r",na = "NA")

# Look at RMA object
CR.Data 

# Automatic population of the Annotation data: 
##### Annotation: mouse4302: corresponds to the GEO entry

# convert exprs to data frame data structure
CR.matrix <- data.frame(exprs(CR.Data))
CR.matrix$affy_id <- rownames(CR.matrix)

head(CR.matrix)   # Data is in log2-scale (required for microarray analysis)

###################################################################################################################
### 2.Retrieve annotation data from BioMart (ensembl) database
###################################################################################################################
library("biomaRt")

# open a connection to biomart
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# see available data in biomart
listMarts(mart, host="www.biomart.org", path="/biomart/martservice", port=80, includeHosts = FALSE, archive=FALSE, verbose = FALSE)
listFilters(mart)
listAttributes(mart)

# choose needed attributes
my.attributes <- c("ensembl_gene_id", "external_gene_name", "affy_mouse430_2")

# get annotations
# load("2020-04-21_Biomart_annotation.RData")  #### result of biomart annotation, save time on webcast day
gene_annot.biomart <- getBM(my.attributes, filters = "affy_mouse430_2", values = CR.matrix$affy_id, mart)
# save(gene_annot.biomart,file = paste0(Sys.Date(),"_Biomart_annotation.RData"))
dim(gene_annot.biomart) #40307    3 

# merge annotation data and Affy data into one data frame
my.annotated.CR.matrix <- merge(gene_annot.biomart, CR.matrix, by.x = 'affy_mouse430_2', by.y='affy_id')

# get informative column names
my.sample.names <- read.table('README_samples.txt', skip = 1, sep = "\t")  # obtained by copy/paste from GEO series page
my.sample.names$celfile <- paste(my.sample.names$V1,".CEL", sep="")

# sanity check
cbind(my.sample.names$celfile,colnames(my.annotated.CR.matrix)[-(1:3)])    # column names correspond !!!
sum(my.sample.names$celfile %in% colnames(my.annotated.CR.matrix)[-(1:3)]) # column names correspond !!!

# rename columns with sample name
colnames(my.annotated.CR.matrix)[-(1:3)] <- my.sample.names$V3

########################
rownames(my.annotated.CR.matrix) <- my.annotated.CR.matrix$external_gene_name
# !!!! does not work : multiple probesets for one gene!!!

# summarize, get max value per gene (or mean, or what makes sense)
my.annotated.CR.matrix.gene <- aggregate(my.annotated.CR.matrix[,-c(1:3)],
                                         by = list(c(my.annotated.CR.matrix$external_gene_name)),
                                         FUN = 'max')
colnames(my.annotated.CR.matrix.gene)[1] <- "Gene_Symbol"
head(my.annotated.CR.matrix.gene)

my.output <- paste(Sys.Date(),"GSE11291_30mths_CR_Data.MaxPerGene_annotated.Data.txt",sep="_")
write.table(my.annotated.CR.matrix.gene, file = my.output, sep = "\t", eol="\r",na = "NA", row.names = F)


###################################################################################################################
### 3. PVclust - hierarchical clustering with bootstrap resampling
###################################################################################################################
library('pvclust')

# call clustering function, enter desired clustering parameters
# For REAL analysis: 1000 or more bootstraps!!!!
result.pvclust <- pvclust(na.omit(my.annotated.CR.matrix[,-(1:3)]), nboot=10) # 10 bootstraps only in the interest of time

# plot results to screen: how are replicates clustering?
plot(result.pvclust)

# plot results to external pdf
pdf(paste(Sys.Date(),"GSE11291_pvclust_postRMA_CR_data.pdf",sep="_"), height = 7, width = 10)
plot(result.pvclust)
dev.off()


###################################################################################################################
### 4. Dimensionality reduction  analysis
###################################################################################################################

# graphical parameters
my.pch.organ <- c(rep(16,20), rep(17,20), rep(18,20))
my.colors.age <- c(rep("coral",5), rep("cyan",5), rep("dodgerblue",5), rep("darkblue",5))

###################################
# a. do MDS analysis
mds.result <- cmdscale(1-cor(my.annotated.CR.matrix.gene[,-(1)],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste(Sys.Date(),"GSE11291_MDS_plot.pdf",sep="_"))
plot(x, y,
     pch = my.pch.organ, col = my.colors.age,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",cex=2)
legend("topright", c("Heart","Gastrocnemius","Neocortex"), col = "grey", pch = c(16,17,18), bty = 'n', pt.cex = 2)
legend("bottomright", c("5m","30m","30m_CR","30m_Resveratrol"), col = c("coral","cyan","dodgerblue","darkblue"), pch = 15, bty = 'n', pt.cex = 2)
dev.off()


###################################
# b. PCA analysis
my.pos.var <- apply(my.annotated.CR.matrix.gene[,-(1)],1,var) > 0
my.pca <- prcomp(t(my.annotated.CR.matrix.gene[my.pos.var,-1]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

pdf(paste(Sys.Date(),"GSE11291_PCA_plot.pdf",sep="_"))
pdf(my.pca.out)
plot(x,y,
     cex=2, pch = my.pch.organ, col = my.colors.age,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5)
legend("topright", c("Heart","Gastrocnemius","Neocortex"), col = "grey", pch = c(16,17,18), bty = 'n', pt.cex = 2)
legend("bottomright", c("5m","30m","30m_CR","30m_Resveratrol"), col = c("coral","cyan","dodgerblue","darkblue"), pch = 15, bty = 'n', pt.cex = 2)
dev.off()


###################################
# c. MDS focusing on heart, "zooming" in to specific conditions of interest
my.heart.CR.data <- my.annotated.CR.matrix.gene[,c(1:21)]
mds.result <- cmdscale(1-cor(my.heart.CR.data[,-(1)],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors <- c(rep("coral",5), rep("cyan",5), rep("dodgerblue",5), rep("darkblue",5))

pdf(paste(Sys.Date(),"GSE11291_Heart_only_MDS_plot.pdf",sep="_"))
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2, pch=16,col=my.colors)
legend("bottomright", c("5m","30m","30m_CR","30m_Resveratrol"), col = c("coral","cyan","dodgerblue","darkblue"), pch = 16, bty = 'n', pt.cex = 2)
dev.off()


###################################
# d. MDS Focusing further on CR impact on 30m heart 
my.heart.CR.data.30m <- my.annotated.CR.matrix.gene[,c(1,7:16)]

# do MDS analysis
mds.result <- cmdscale(1-cor(my.heart.CR.data.30m[,-(1)],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors <- c(rep("cyan",5), rep("dodgerblue",5))

pdf(paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_only_MDS_plot.pdf",sep="_"))
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2, pch=16,col=my.colors)
legend("topright",c("30m_CTL","30m_CR"),col=c("cyan","dodgerblue"),pch=16,bty='n',pt.cex=2)
dev.off()


###################################################################################################################
#### 5. Differential gene expression analysis with Limma
###################################################################################################################
library(limma)

my.diet <- factor(c(rep("CTL",5),rep("CR",5)), levels=c("CTL","CR"))
rownames(my.heart.CR.data.30m) <- my.heart.CR.data.30m$Gene_Symbol

# build model matrix
my.model.heart <- model.matrix( ~ my.diet)

# fit the limma model (generalization of a linear model for microarray data)
my.fit.heart <- lmFit(my.heart.CR.data.30m[,-c(1)], my.model.heart)
my.fit.heart <- eBayes(my.fit.heart)

my.sig.heart <- topTable(my.fit.heart, coef = "my.dietCR", p.value = 1 , number = Inf)
save(my.sig.heart, file = paste(Sys.Date(),"GSE11291_Heart_CRvsCTL.RData",sep="_"))

# Diagnostic plots
volcanoplot(my.fit.heart, coef = "my.dietCR")
plotMA(my.fit.heart, coef = "my.dietCR")

# manual volcano plot with significant genes highlighted
pdf(paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_limma_volcano_Plot.pdf",sep="_"))
plot(my.sig.heart$logFC, -log10(my.sig.heart$adj.P.Val), 
     pch = 16, cex = 0.5,
     xlab = "log2(Fold Change CR/CTL)",
     ylab = "-log10(FDR)",
     xlim = c(-4,4), ylim = c(0,10))
points(my.sig.heart$logFC[my.sig.heart$adj.P.Val < 0.05], -log10(my.sig.heart$adj.P.Val[my.sig.heart$adj.P.Val < 0.05]),
       col = "red",pch = 16, cex = 0.5)
legend("topleft",c("FDR < 5%"),col=c("red"),pch=16,bty='n',pt.cex=0.5)
dev.off()


pdf(paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_limma_volcano_Plot_smoothScatter.pdf",sep="_"))
smoothScatter(my.sig.heart$logFC, -log10(my.sig.heart$adj.P.Val), 
              xlab = "log2(Fold Change CR/CTL)",
              ylab = "-log10(FDR)")
points(my.sig.heart$logFC[my.sig.heart$adj.P.Val < 0.05], -log10(my.sig.heart$adj.P.Val[my.sig.heart$adj.P.Val < 0.05]),
       col = "red",pch = 16, cex = 0.5,
       xlim = c(-4,4), ylim = c(0,10))
legend("topleft",c("FDR < 5%"),col=c("red"),pch=16,bty='n',pt.cex=0.5)
dev.off()

my.output.heart <- paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_limma_All_results.txt",sep="_")
write.table(my.sig.heart, file = my.output.heart, sep = "\t", quote = F)


###################################################################################################################
#### 6. Heatmap of differentially expressed genes
##################################################################################################################

### get the heatmap of aging changes at FDR5
my.sig.heart <- my.sig.heart[!is.na(my.sig.heart$adj.P.Val),] ## need to exclude NA

genes.heart.CR <- rownames(my.sig.heart)[my.sig.heart$adj.P.Val < 0.05] # FDR < 5%

my.num.heart.CR <- length(genes.heart.CR)
my.num.heart.CR

# heatmap drawing 
library('pheatmap')
my.heatmap.out <- paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_Heatmap_significant_genes.pdf",sep="_")

pdf(my.heatmap.out, height = 20, width = 10, onefile = F)
my.heatmap.title <- paste("Heart CR significant (FDR < 5%), ",my.num.heart.CR, " genes",sep="")
pheatmap(my.heart.CR.data.30m[genes.heart.CR,-1],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 30)
dev.off()



###################################################################################################################
#### 7. get DE genes from CR in gastrocnemius for comparison (repeat of above on new data)
##################################################################################################################

# extract gastrocnemius data
my.gastroc.CR.data.30m <- my.annotated.CR.matrix.gene[,c(1,27:36)]

# do MDS analysis on gastrocnemius
mds.result <- cmdscale(1-cor(my.gastroc.CR.data.30m[,-(1)],method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste(Sys.Date(),"GSE11291_gastrocnemius_CRvsCTL_only_MDS_plot.pdf",sep="_"))
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2, pch=16,col=my.colors)
legend("topleft",c("30m_CTL","30m_CR"),col=c("cyan","dodgerblue"),pch=16,bty='n',pt.cex=2)
dev.off()

# run limma analysis
my.diet <- factor(c(rep("CTL",5),rep("CR",5)), levels=c("CTL","CR"))    # same as above, repeated for readability
rownames(my.gastroc.CR.data.30m) <- my.gastroc.CR.data.30m$Gene_Symbol

# build model matrix
my.model.gastroc <- model.matrix( ~ my.diet)

# fit the limma model (generalization of a linear model for microarray data)
my.fit.gastroc  <- lmFit(my.gastroc.CR.data.30m[,-c(1)], my.model.gastroc)
my.fit.gastroc  <- eBayes(my.fit.gastroc)

my.sig.gastroc  <- topTable(my.fit.gastroc , coef = "my.dietCR", p.value = 1 , number = Inf)
save(my.sig.gastroc , file = paste(Sys.Date(),"GSE11291_Gastrocnemius_CRvsCTL.RData",sep="_"))

# manual volcano plot with significant genes highlighted
pdf(paste(Sys.Date(),"GSE11291_Gastrocnemius_CRvsCTL_limma_volcano_Plot.pdf",sep="_"))
plot(my.sig.gastroc$logFC, -log10(my.sig.gastroc$adj.P.Val), 
     pch = 16, cex = 0.5,
     xlab = "log2(Fold Change CR/CTL)",
     ylab = "-log10(FDR)",
     xlim = c(-4,4), ylim = c(0,10))
points(my.sig.gastroc$logFC[my.sig.gastroc$adj.P.Val < 0.05], -log10(my.sig.gastroc$adj.P.Val[my.sig.gastroc$adj.P.Val < 0.05]),
       col = "red",pch = 16, cex = 0.5)
legend("topleft",c("FDR < 5%"),col=c("red"),pch=16,bty='n',pt.cex=0.5)
dev.off()

my.output <- paste(Sys.Date(),"GSE11291_Gastrocnemius_CRvsCTL_limma_All_results.txt",sep="_")
write.table(my.sig.gastroc, file = my.output, sep = "\t", quote = F)


### get the heatmap of aging changes at FDR5
my.sig.gastroc <- my.sig.gastroc[!is.na(my.sig.gastroc$adj.P.Val),] ## need to exclude NA
genes.gastroc.CR <- rownames(my.sig.gastroc)[my.sig.gastroc$adj.P.Val < 0.05] # FDR < 5%
my.num.gastroc.CR <- length(genes.gastroc.CR)
my.num.gastroc.CR

# heatmap drawing 
my.heatmap.out <- paste(Sys.Date(),"GSE11291_Heart_CRvsCTL_Heatmap_significant_genes.pdf",sep="_")
pdf(my.heatmap.out, height = 20, width = 10, onefile = F)
my.heatmap.title <- paste("Heart CR significant (FDR < 5%), ",my.num.gastroc.CR, " genes",sep="")
pheatmap(my.gastroc.CR.data.30m[genes.gastroc.CR,-1],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 30)
dev.off()


###################################################################################################################
#### 8. compare DE genes under CR in heart vs. gastrocnemius
##################################################################################################################
library(bitops)

my.heart.up   <- rownames(my.sig.heart)[bitAnd(my.sig.heart$logFC > 0, my.sig.heart$adj.P.Val < 0.05)>0]
my.gastroc.up <- rownames(my.sig.gastroc)[bitAnd(my.sig.gastroc$logFC > 0, my.sig.gastroc$adj.P.Val < 0.05)>0]

# get Venn Diagram
library(Vennerable)

my.up.genes <- list("Heart"            = my.heart.up,
                    "Gastrocnemius"    = my.gastroc.up)
my.Venn <- Venn(my.up.genes)

my.venn.out <- paste(Sys.Date(),"GSE11291_Heart_vs_Gastroc_CR_Venn.pdf",sep="_")
pdf(my.venn.out, height = 20, width = 10, onefile = F)
plot(my.Venn, doWeights=T, show=list(Faces=FALSE))
dev.off()

# is the overlap significant?

# get the "universe/background"
my.background <- intersect(rownames(my.sig.heart),rownames(my.sig.gastroc))

my.contingency.mat <- matrix(c(618,3086,2446,length(my.background)-618-3086-2446),2,2)

# perform Fisher exact/hypergeometric test
fisher.test(my.contingency.mat, alternative = "greater")
# p-value = 2.3e-10
# overlap is greater than expected by chance





