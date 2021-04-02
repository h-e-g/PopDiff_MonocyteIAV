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
############ with and without stimulation, including T0                            ###############################
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

#isolate the genes with average expression > 0.1 in at least one subset plus condition at T0-T8
keepVec <- avgCondSub[do.call(pmax, c(avgCondSub[2:19], list(na.rm=TRUE))) > 0.1, 'gene']
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

#one master analysis
twoSubs <- findMarkers(mat[keepVec,], keepMeta$sxCD16, test.type=c("wilcox"), block=paste(keepMeta$COND, keepMeta$TP, keepMeta$ID1, sep="_"))
twoSubs[[1]]$gene <- rownames(twoSubs[[1]])
twoSubsT <- findMarkers(mat[keepVec,], keepMeta$sxCD16, test.type=c("t"), block=paste(keepMeta$COND, keepMeta$TP, keepMeta$ID1, sep="_"))
twoSubs[[1]]$logFC <- twoSubsT[[1]][match(twoSubs[[1]]$gene, rownames(twoSubsT[[1]])), 4]
ts <- data.frame(twoSubs[[1]][,c(5, 6, 4, 2, 3)])
names(ts) <- c("gene", 
               "CD16neg_CD16pos_logFC",
               "CD16neg_CD16pos_AUC",
               "CD16neg_CD16pos_pval",
               "CD16neg_CD16pos_fdr")

#What genes are stable & interesting? We again use a logFC of 0.2 and FDR=1
ts$color <- ifelse(ts$CD16neg_CD16pos_fdr >= 0.01 | abs(ts$CD16neg_CD16pos_logFC) < 0.2, 'nonsig',
                   ifelse(ts$CD16neg_CD16pos_logFC > 0, "CD16-", "CD16+"))

#Add ensembl ids for genes
ts$ens <- GA[match(ts$gene, GA$geneAlt), 'stable']
which( is.na( ts$ens ) )

#Combine
master <- right_join(avgCondSub, t0, by="gene")
master <- right_join(master, ni, by="gene")
master <- left_join(master, iav, by="gene")
master <- left_join(master, ts, by="gene")

#Stats
cor.test(master$ni_CD16neg_CD16pos_logFC, master$iav_CD16neg_CD16pos_logFC)

logFC.ni <- data.frame(T0=master$T0_NI_CD16pos-master$T0_NI_CD16neg,
                       T2=master$T2_NI_CD16pos-master$T2_NI_CD16neg, 
                       T4=master$T4_NI_CD16pos-master$T4_NI_CD16neg, 
                       T6=master$T6_NI_CD16pos-master$T6_NI_CD16neg, 
                       T8=master$T8_NI_CD16pos-master$T8_NI_CD16neg)
cor(logFC.ni)

logFC.iav <- data.frame(T0=master$T0_NI_CD16pos-master$T0_NI_CD16neg,
                       T2=master$T2_IAV_CD16pos-master$T2_IAV_CD16neg, 
                       T4=master$T4_IAV_CD16pos-master$T4_IAV_CD16neg, 
                       T6=master$T6_IAV_CD16pos-master$T6_IAV_CD16neg, 
                       T8=master$T8_IAV_CD16pos-master$T8_IAV_CD16neg)
cor(logFC.iav)

#save the data frame
write.table(master, file=paste0(PopDiff, "results/T0-T8_twoSubsets_5681genes.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#volcano plot
ggplot(ts) +
  geom_point(aes(x=CD16neg_CD16pos_logFC, y=-log10(CD16neg_CD16pos_fdr), color=color), alpha=0.5, size=0.5) +
  theme_light() +
  scale_color_manual(values=c("#666666", "#000000", "#CCCCCC")) +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(),
        legend.title=element_blank(),
        legend.position="none",
        text = element_text(size=8),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  xlab("logFC") +
  ylab("-log10(FDR)")

#Define gene lists for enrichment analyses
CD16neg <- ts[ts$color=="CD16-", 'ens']
CD16pos <- ts[ts$color=="CD16+", 'ens']

source(paste0(PopDiff, "/scripts/GOSeq_hg38.R"))

background = ts$ens

enrLists = list(CD16neg, CD16pos)
names(enrLists) <- c("CD16neg", "CD16pos")

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

write.table(enr[,-10], paste0(PopDiff, "results/CD16pos_CD16neg_byTPCOND_GO.enrichments_5681backgroundgenes.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#write.table(enrLists[["CD16neg"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16neg_IAV_308genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
#write.table(enrLists[["CD16pos"]], file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/CD16pos_IAV_400genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))
#write.table(background, file=paste0(PopDiff, "results/gene_lists/T2-T8_IAV_NI_Subsets/background_5368genes.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))

##################################################################################################################
############ End #################################################################################################
##################################################################################################################
