#Setup - Make sure the following programs are installed and loaded
library(scran)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(reshape2)
library(Matrix)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#Define the working directory - Alter to configure with your system
PopDiff="/Users/moneill/Desktop/PopDiff_MonocyteIAV/"

#Load the cleaned and normalized single-cell RNA-seq data
keepMeta=readRDS(paste0(PopDiff, "data/meta_dataframe_88559monocytes.rds"))
mat=readRDS(paste0(PopDiff, "/data/normalized_count_matrix_88559monocytes_22593hostgenes.rds"))
GA=readRDS(paste0(PopDiff, "data/gene_annotation_22593hostgenes.rds"))

##################################################################################################################
############ Calculate differentially expressed genes between CD16-/CD16+ subsets  ###############################
############ with and without stimulation                                          ###############################
##################################################################################################################
identical(colnames(mat), keepMeta$Barcode)

# calculate average gene expression per CD16+/- subset per condition
tmpVar <- paste(keepMeta$TP_COND, keepMeta$sxCD16, sep="_")
tmpList = list()
for (i in levels(as.factor(tmpVar))){
  avg <- unname(rowMeans(mat[, tmpVar==i]))
  tmpList[[i]] = avg
}
avgCondSub <- data.frame(tmpList)
avgCondSub$gene <- rownames(mat)
avgCondSub <- avgCondSub[,c(19,1:18)]
colnames(avgCondSub) <- gsub("X", "T", colnames(avgCondSub))
colnames(avgCondSub) <- gsub("FALSE", "CD16neg", colnames(avgCondSub))
colnames(avgCondSub) <- gsub("TRUE", "CD16pos", colnames(avgCondSub))

#isolate the genes with average expression > 0.1 in at least one subset plus condition at T2-T8
keepVec <- avgCondSub[do.call(pmax, c(avgCondSub[4:19], list(na.rm=TRUE))) > 0.1, 'gene']
length(keepVec)

#T2-T8 IAV
subsetVar <- keepMeta$TP!=0 & keepMeta$COND == "IAV"

iavW <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("wilcox"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
iavW[[1]]$gene <- rownames(iavW[[1]])
iavT <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("t"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
iavW[[1]]$logFC <- iavT[[1]][match(iavW[[1]]$gene, rownames(iavT[[1]])), 4]
iav <- data.frame(iavW[[1]][,c(5, 6, 4, 2, 3)])
names(iav) <- c("gene", 
                "iav_CD16neg_CD16pos_logFC",
                "iav_CD16neg_CD16pos_AUC",
                "iav_CD16neg_CD16pos_pval",
                "iav_CD16neg_CD16pos_fdr")

#T2-T8 NI
subsetVar <- keepMeta$TP!=0 & keepMeta$COND == "NI"

niW <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("wilcox"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
niW[[1]]$gene <- rownames(niW[[1]])
niT <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("t"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
niW[[1]]$logFC <- niT[[1]][match(niW[[1]]$gene, rownames(niT[[1]])), 4]
ni <- data.frame(niW[[1]][,c(5, 6, 4, 2, 3)])
names(ni) <- c("gene", 
               "ni_CD16neg_CD16pos_logFC",
               "ni_CD16neg_CD16pos_AUC",
               "ni_CD16neg_CD16pos_pval",
               "ni_CD16neg_CD16pos_fdr")

#T0 (for Comparison)
subsetVar <- keepMeta$TP==0

t0W <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("wilcox"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
t0W[[1]]$gene <- rownames(t0W[[1]])
t0T <- findMarkers(mat[keepVec,subsetVar], keepMeta[subsetVar, 'sxCD16'], test.type=c("t"), block=paste(keepMeta[subsetVar, 'TP'], keepMeta[subsetVar, 'ID1'], sep="_"))
t0W[[1]]$logFC <- t0T[[1]][match(t0W[[1]]$gene, rownames(t0T[[1]])), 4]
t0 <- data.frame(t0W[[1]][,c(5, 6, 4, 2, 3)])
names(t0) <- c("gene", 
               "t0_CD16neg_CD16pos_logFC",
               "t0_CD16neg_CD16pos_AUC",
               "t0_CD16neg_CD16pos_pval",
               "t0_CD16neg_CD16pos_fdr")

#Combine
t2t8.master <- right_join(avgCondSub, ni, by="gene")
t2t8.master <- left_join(t2t8.master, iav, by="gene")

#What genes are interesting? We again use a logFC of 0.2 and FDR=1
t2t8.master$IAV <- ifelse(t2t8.master$iav_CD16neg_CD16pos_fdr >= 0.01 | abs(t2t8.master$iav_CD16neg_CD16pos_logFC) < 0.2, 'nonsig',
                        ifelse(t2t8.master$iav_CD16neg_CD16pos_logFC > 0, "classical", "nonclassical"))
t2t8.master$NI <- ifelse(t2t8.master$ni_CD16neg_CD16pos_fdr >= 0.01 | abs(t2t8.master$ni_CD16neg_CD16pos_logFC) < 0.2, 'nonsig',
                          ifelse(t2t8.master$ni_CD16neg_CD16pos_logFC > 0, "classical", "nonclassical"))

#Add ensembl ids for genes
t2t8.master$ens <- GA[match(t2t8.master$gene, GA$geneAlt), 'stable']
which( is.na( t2t8.master$ens ) )

#Add the T0 DE analysis based on CD16-/CD16+ subsets
t2t8.master <- left_join(t2t8.master, t0, by="gene")

#save the data frame
write.table(t2t8.master, file=paste0(PopDiff, "results/postT0_CD16neg_CD16pos_byCond_5368genes.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#Define gene lists for enrichment analyses
CD16neg_IAV <- t2t8.master[t2t8.master$IAV=="classical", 'ens']
CD16pos_IAV <- t2t8.master[t2t8.master$IAV=="nonclassical", 'ens']
CD16neg_NI <- t2t8.master[t2t8.master$NI=="classical", 'ens']
CD16pos_NI <- t2t8.master[t2t8.master$NI=="nonclassical", 'ens']

source(paste0(PopDiff, "/scripts/GOSeq_hg38.R"))

background = t2t8.master$ens

enrLists = list(CD16neg_IAV, CD16pos_IAV, CD16neg_NI, CD16pos_NI)
names(enrLists) <- c("CD16neg_IAV", "CD16pos_IAV", "CD16neg_NI", "CD16pos_NI")

paste("Performing enrichment")
subGO=list()
for (y in 1:length(enrLists)){
  t=unique(enrLists[[y]])
  goRES <- GOSeq(t, background, addCI=TRUE)
  if (length(goRES$category) > 0) {
    goRES$genes <- NA
    tmpGA <- dplyr::filter(GA, stable %in% t)
    for (i in 1:length(goRES$category)) {
      category=goRES$category[i]
      genes <- tmpGA[grepl(category, tmpGA$GO),'Symbol']
      genes <- paste(genes, collapse=",")
      goRES[i, 'genes'] = genes
    }
  }
  nm=names(enrLists)[y]
  subGO[[nm]]=goRES
}

enr <- rbindlist(subGO, idcol=T, fill=T)

write.table(enr[,-10], paste0(PopDiff, "results/postT0_CD16pos_CD16neg_byCOND_GO.enrichments_5368backgroundgenes.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

write.table(enrLists[["CD16neg_IAV"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16neg_IAV_308genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(enrLists[["CD16pos_IAV"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16pos_IAV_400genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(enrLists[["CD16neg_NI"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16neg_NI_280genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(enrLists[["CD16pos_NI"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16pos_NI_402genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
write.table(background, file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/background_5368genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))

##################################################################################################################
############ End #################################################################################################
##################################################################################################################
