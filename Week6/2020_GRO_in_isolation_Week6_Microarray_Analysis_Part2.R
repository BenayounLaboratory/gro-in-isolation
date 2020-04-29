##########################################
############ GRO in isolation ############
##########################################
# Copyright 2020, Bérénice A. Benayoun
# 2020-04-29

#######        Week 6 code         #######

### Example of analysis of a microarray experiment  (Part 2) ###
### Functional annotation enrichment analysis

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
setwd('/Users/berenice/Dropbox/USC_Work_Folder/Teaching/2019-2020/GRO-in-isolation/Github/Week6')
options(stringsAsFactors = F)

# read in key files from last week's session
my.limma.gastroc <- read.table('2020-04-22_GSE11291_Gastrocnemius_CRvsCTL_limma_All_results.txt', header = T, sep = "\t")
my.norm.data     <- read.table('2020-04-22_GSE11291_30mths_CR_Data.MaxPerGene_annotated.Data.txt', header = T, sep = "\t")


###################################################################################################################
#### 9. Functional enrichment analysis by Over Representation Analysis (ORA) - clusterProfiler
###################################################################################################################
library(org.Mm.eg.db)
library(clusterProfiler)

# extract genes up/downregulated with FDR < 5% from limma results
gastroc.CR.up <- intersect(rownames(my.limma.gastroc)[my.limma.gastroc$adj.P.Val < 0.05], 
                           rownames(my.limma.gastroc)[my.limma.gastroc$logFC > 0])
gastroc.CR.dwn <- intersect(rownames(my.limma.gastroc)[my.limma.gastroc$adj.P.Val < 0.05], 
                            rownames(my.limma.gastroc)[my.limma.gastroc$logFC < 0])
gastroc.background <- rownames(my.limma.gastroc)

# Convert to ENTREZ ID (necessary for on board KEGG analysis)
# keytypes(org.Mm.eg.db)
entrezID.gastroc.CR.up       <- bitr(gastroc.CR.up     , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.gastroc.CR.dwn      <- bitr(gastroc.CR.dwn    , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.gastroc.background  <- bitr(gastroc.background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


##########################################################
########    A. GO BP over-representation test     ########
ego.bp.up <- enrichGO(gene          = entrezID.gastroc.CR.up$ENTREZID,
                      universe      = entrezID.gastroc.background$ENTREZID,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      minGSSize     = 100,   # just for speed purposes here!
                      maxGSSize     = 1000,   # just for speed purposes here!
                      readable      = TRUE)

# ego.bp.dwn <- enrichGO(gene          = entrezID.gastroc.CR.dwn$ENTREZID,
#                        universe      = entrezID.gastroc.background$ENTREZID,
#                        OrgDb         = org.Mm.eg.db,
#                        keyType       = 'ENTREZID',
#                        ont           = "BP",
#                        pAdjustMethod = "BH",
#                        pvalueCutoff  = 0.01,
#                        qvalueCutoff  = 0.05,
#                        readable      = TRUE)

# write results to file
write.table(ego.bp.up@result,  file = paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
# write.table(ego.bp.dwn@result, file = paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")

my.GOBP.res <- ego.bp.up@result

# make some plots
pdf(paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_dotplot_FDR5_UP.pdf", sep = "_"))
dotplot(ego.bp.up,  x = "Count", title = "Upregulated GO BP Terms")
dev.off()

pdf(paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_GeneConcept_Network_FDR5_UP.pdf", sep = "_"))
cnetplot(ego.bp.up,  categorySize="pvalue", node_label = 'category')
dev.off()

pdf(paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_Enrichment_Map_FDR5_UP.pdf", sep = "_"))
emapplot(ego.bp.up)
dev.off()

##########################################################
########     B. KEGG over-representation test     ########
kk.up  <- enrichKEGG(gene          = entrezID.gastroc.CR.up$ENTREZID,
                     universe      = entrezID.gastroc.background$ENTREZID,
                     keyType       = 'ncbi-geneid',
                     organism      = 'mmu',
                     pvalueCutoff  = 0.05)

# kk.dwn <- enrichKEGG(gene          = entrezID.gastroc.CR.dwn$ENTREZID,
#                      universe      = entrezID.gastroc.background$ENTREZID,
#                      keyType       = 'ncbi-geneid',
#                      organism      = 'mmu',
#                      pvalueCutoff  = 0.05)

# write results to file
write.table(kk.up@result,  file = paste(Sys.Date(),"Gastrocnemius_CR_KEGG_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
# write.table(kk.dwn@result, file = paste(Sys.Date(),"Gastrocnemius_CR_KEGG_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")

my.kegg.res <- kk.up@result


# make some plots
pdf(paste(Sys.Date(),"Gastrocnemius_CR_KEGG_dotplot_FDR5_UP.pdf", sep = "_"))
dotplot(kk.up,  title = "Upregulated KEGG Pathways")
dev.off()

pdf(paste(Sys.Date(),"Gastrocnemius_CR_KEGG_GeneConcept_Network_FDR5_UP.pdf", sep = "_"))
cnetplot(kk.up,  categorySize="pvalue", node_label = 'category')
dev.off()

pdf(paste(Sys.Date(),"Gastrocnemius_CR_GO_BP_Enrichment_Map_FDR5_UP.pdf", sep = "_"))
emapplot(kk.up)
dev.off()

###################################################################################################################
#### 10. Functional enrichment analysis by Gene Set Enrichment Analysis (GSEA) - clusterProfiler
###################################################################################################################
library(org.Mm.eg.db)
library(clusterProfiler)

# Prepare GeneList using Limma t-statistic to rank genes
my.limma.gastroc$symbol    <- rownames(my.limma.gastroc)
my.limma.gastroc.v2        <- merge(data.frame(my.limma.gastroc), entrezID.gastroc.background, by.x = "symbol", by.y = "SYMBOL")
gastroc.CR.geneList        <- my.limma.gastroc.v2$t
names(gastroc.CR.geneList) <- my.limma.gastroc.v2$ENTREZID
gastroc.CR.geneList        <- sort(gastroc.CR.geneList, decreasing = TRUE)


##########################################################
######## A. GO Gene Set Enrichment Analysis
go.bp.gsea <- gseGO(geneList     = gastroc.CR.geneList,
                    OrgDb        = org.Mm.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "BP",
                    nPerm        = 1000,
                    minGSSize    = 100,   # in the interest of time
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)

# write results to file
write.table(go.bp.gsea@result, file = paste(Sys.Date(),"Gastrocnemius_CR_GOBP_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
my.go.bp.gsea.res <- go.bp.gsea@result

# GSEA plot of interesting pathways
gseaplot(go.bp.gsea, geneSetID = "GO:0016570", title = "GO:0016570 histone modification")                      # downregulated example
gseaplot(go.bp.gsea, geneSetID = "GO:0006281", title = "GO:0006281 DNA repair")
gseaplot(go.bp.gsea, geneSetID = "GO:0007606", title = "GO:0007606 sensory perception of chemical stimulus")   # upregulated example


##########################################################
######## B. KEGG Gene Set Enrichment Analysis
kegg.gsea <- gseKEGG(geneList     = gastroc.CR.geneList,
                     organism     = 'mmu',
                     keyType      = 'ncbi-geneid',
                     nPerm        = 1000,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)

# write results to file
write.table(kegg.gsea@result, file = paste(Sys.Date(),"Gastrocnemius_CR_KEGG_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
my.kegg.gsea.res <- kegg.gsea@result

# GSEA plot of interesting pathways
gseaplot(kegg.gsea, geneSetID = "mmu04151", title = "PI3K-Akt signaling pathway")                      # downregulated example
gseaplot(kegg.gsea, geneSetID = "mmu00591", title = "Linoleic acid metabolism")                    # upregulated example


###################################################################################################################
#### 10. A completely different approach : Gene Set Variation Analysis (GSVA)
###################################################################################################################
library(limma)
library(GSVA)
library(GSEABase)

#################################################################
######### Functions to make GMT files "mouse-friendly"  ######### 
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

convert_to_mouse <- function(my.coll) {
  for (i in 1:length(my.coll@.Data)) {
    my.coll@.Data[[i]]@geneIds <- firstup( tolower( my.coll@.Data[[i]]@geneIds) )
  }
  return(my.coll)
}
##################################################################

# read in geneset annotation file in gmt format (Here KEGG data from MSigDB)
my.c2.kegg  <- convert_to_mouse(getGmt("c2.cp.kegg.v7.0.symbols.gmt", geneIdType = SymbolIdentifier(), collectionType =  KEGGCollection() ))

# see MSigDB and Harmonizome webiste to get these files
# MSigDB      : https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Harmonizome : http://amp.pharm.mssm.edu/Harmonizome/download
# You can also make custom files!!!!

# extract gastrocnemius CR data
rownames(my.norm.data) <- my.norm.data$Gene_Symbol
my.gastroc.CR.data.30m <- my.norm.data[,c(27:36)]

# build model matrix
my.diet <- factor(c(rep("CTL",5),rep("CR",5)), levels=c("CTL","CR"))    # same as above, repeated for readability
my.model.gastroc <- model.matrix( ~ my.diet)

# run GSVA to get geneset activity in each sample
my.exp.coll <- gsva(as.matrix(my.gastroc.CR.data.30m), # needs matrix input
                    my.c2.kegg,
                    method="gsva",
                    kcdf="Gaussian",
                    min.sz=10,
                    max.sz=500)
head(my.exp.coll)

# run limma fit on GSVA matrix
my.fit.coll <- lmFit(my.exp.coll, design = my.model.gastroc)
my.fit.coll <- eBayes(my.fit.coll)

my.all.coll <- topTable(my.fit.coll, coef="my.dietCR", number= Inf )
my.DE.coll  <- topTable(my.fit.coll, coef="my.dietCR", number= Inf, p.value = 0.05, adjust="BH") # 95 significant pathways

write.table(my.DE.coll, file = paste0(Sys.Date(),"_KEGG_GSVA_Analysis_GeneSets_FDR5.txt"), sep = "\t", quote = F)

my.outprefix <- paste(Sys.Date(),"KEGG_GSVA_Analysis_GeneSets_FDR5", sep = "_")

# output heatmap of top 15 significant pathways
library('pheatmap')
pheatmap(my.exp.coll[rownames(my.DE.coll)[1:10],], 
         scale = 'row',
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         border_color = NA,
         show_rownames = T,
         cellwidth = 20,
         cellheight = 10)



###################################################################################################################
#### 11. Functional enrichment analysis: outside R
###################################################################################################################

##### Online tools
# STRING
# EnrichR
# DAVID

# Motif analysis in promoters?


