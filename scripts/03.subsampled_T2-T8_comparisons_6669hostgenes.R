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
library(dynamicTreeCut)
library(ComplexHeatmap)
library(inlmisc)

#Define the working directory - Alter to configure with your system
PopDiff="/Users/moneill/Desktop/PopDiff_MonocyteIAV/"

#Load the cleaned and normalized single-cell RNA-seq data
keepMeta=readRDS(paste0(PopDiff, "data/meta_dataframe_88559monocytes.rds"))
mat=readRDS(paste0(PopDiff, "/data/normalized_count_matrix_88559monocytes_22593hostgenes.rds"))
GA=readRDS(paste0(PopDiff, "data/gene_annotation_22593hostgenes.rds"))

##################################################################################################################
############ Calculate differentially expressed genes on subsampled dataset  #####################################
##################################################################################################################
#Treat SoupX CD16 calls as a factor
keepMeta$sxCD16 <- as.factor(keepMeta$sxCD16)

#Check that meta and matrix have same order
identical(keepMeta$Barcode, colnames(mat))

#Randomly subsample 100 cells from each TP+status
set.seed(1000) #for reproducibility
dwnsmp <- do.call( rbind, lapply( split(keepMeta, keepMeta$c30) ,
                                  function(df) {
                                    idweights = data.frame(1/(table(df$ID1)/sum(table(df$ID1))))
                                    probs <- idweights[match(df$ID1, idweights$Var1), 'Freq'] #increase prob of sampling to donors with fewer cells
                                    subSamp <- df[sample(nrow(df), 100, prob=probs), ]
                                    return(subSamp)
                                  }))

dwnsmp$cluster <- gsub("low", "infected", dwnsmp$cluster)
dwnsmp$cluster <- gsub("high", "infected", dwnsmp$cluster)
table(dwnsmp$cluster)

# Randomly subsample 100 cells at T4 for each infected cluster
t4Inf <- keepMeta[keepMeta$TP == 4 & keepMeta$COND == "IAV" & keepMeta$sxFlu, ]
t4Inf$c10 <- as.factor(as.character(t4Inf$c10))  

set.seed(1000)
t4Inf <- do.call( rbind, lapply( split(t4Inf, t4Inf$c10) ,
                                 function(df) {
                                   idweights = data.frame(1/(table(df$ID1)/sum(table(df$ID1))))
                                   probs <- idweights[match(df$ID1, idweights$Var1), 'Freq']
                                   subSamp <- df[sample(nrow(df), 100, prob=probs), ]
                                   return(subSamp)
                                 }))

table(t4Inf$cluster)
t4df <- rbind(t4Inf, filter(dwnsmp, TP==4 & cluster=="bystander")) #combine with 200 bystander cells
table(t4df$cluster)

#Find the mean for each cluster
dwn <- mat[,dwnsmp$Barcode]
t4 <- mat[,t4Inf$Barcode]
t4v2 <- mat[,t4df$Barcode]

tmpList = list()
for (i in levels(as.factor(as.character(dwnsmp$c30)))){
  avg <- unname(rowMeans(dwn[,dwnsmp$c30==i]))
  tmpList[[i]] = avg
}
tmpDF <- data.frame(tmpList)
names(tmpDF) <- gsub("X", "T", names(tmpDF))
names(tmpDF) <- gsub("inf", "infected", names(tmpDF))
names(tmpDF) <- gsub("TRUE", "CD16pos", names(tmpDF))
names(tmpDF) <- gsub("FALSE", "CD16neg", names(tmpDF))

master <- cbind(gene=rownames(dwn), tmpDF)

#add mean for T4 high vs low
tmpList = list()
for (i in levels(as.factor(as.character(t4Inf$c34)))){
  avg <- unname(rowMeans(t4[,t4Inf$c34==i]))
  tmpList[[i]] = avg
}
tmpDF <- data.frame(tmpList)
names(tmpDF) <- gsub("X", "T", names(tmpDF))
names(tmpDF) <- gsub("TRUE", "CD16pos", names(tmpDF))
names(tmpDF) <- gsub("FALSE", "CD16neg", names(tmpDF))

master <- cbind(master, tmpDF)

#Define function to run the cond/fate comparisons
fate_comp <- function(state1, state2) {
  
  #CD16neg
  print("CD16neg subset")
  subsetVar <- (dwnsmp$cluster == state1 | dwnsmp$cluster == state2) & dwnsmp$sxCD16=="FALSE"
  modDS <- dwnsmp[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  CD16neg <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"), block=modDS$TP)
  CD16neg[[1]]$gene <- rownames(CD16neg[[1]])
  logFC <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("t"), block=modDS$TP)
  CD16neg[[1]]$logFC <- logFC[[1]][match(CD16neg[[1]]$gene, rownames(logFC[[1]])), 4]
  CD16neg <- data.frame(CD16neg[[1]][,c(5, 6, 4, 2, 3)])
  names(CD16neg) <- c("gene",
                      paste("CD16neg", state1, state2, "logFC", sep="_"),
                      paste("CD16neg", state1, state2, "AUC", sep="_"),
                      paste("CD16neg", state1, state2, "pval", sep="_"), 
                      paste("CD16neg", state1, state2, "fdr", sep="_"))
  head(CD16neg)
  
  #CD16pos
  print("CD16pos subset")
  subsetVar <- (dwnsmp$cluster == state1 | dwnsmp$cluster == state2) & dwnsmp$sxCD16=="TRUE"
  modDS <- dwnsmp[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  CD16pos <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"), block=modDS$TP)
  CD16pos[[1]]$gene <- rownames(CD16pos[[1]])
  logFC <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("t"), block=modDS$TP)
  CD16pos[[1]]$logFC <- logFC[[1]][match(CD16pos[[1]]$gene, rownames(logFC[[1]])), 4]
  CD16pos <- data.frame(CD16pos[[1]][,c(5, 6, 4, 2, 3)])
  names(CD16pos) <- c("gene",
                      paste("CD16pos", state1, state2, "logFC", sep="_"),
                      paste("CD16pos", state1, state2, "AUC", sep="_"),
                      paste("CD16pos", state1, state2, "pval", sep="_"), 
                      paste("CD16pos", state1, state2, "fdr", sep="_"))
  head(CD16pos)
  
  
  #Subset differences
  print("Testing for differences between subsets")
  subsetVar <- (dwnsmp$cluster == state1 | dwnsmp$cluster == state2)
  modDS <- dwnsmp[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  dataList=list(apply(dwn[,subsetVar], 1, function(y){
    mod <- lm(y ~ modDS$sxCD16 + modDS$cluster + modDS$sxCD16:modDS$cluster + as.factor(modDS$TP))
    c(summary(mod)$coefficients[7,4])}))
  df <- data.frame(dataList)
  df$gene <- rownames(df)
  df <- df[,c(2,1)]
  names(df) <- c("gene", "pval")
  df$fdr <- p.adjust(df$pval, method='fdr')
  names(df) <- c("gene", paste("subsetDiff", state1, state2, "pval", sep="_"), paste("subsetDiff", state1, state2, "fdr", sep="_"))
  head(df)
  
  ofInt <- full_join(CD16pos, CD16neg, by="gene")
  ofInt <- full_join(ofInt, df, by="gene")
  
  return(ofInt)
}

#Run function on each comparison
un.by <- fate_comp("unexposed", "bystander")
un.inf <- fate_comp("unexposed", "infected")
by.inf <- fate_comp("bystander", "infected")

#Define function to run the cond/fate comparisons, combining subsets *************** (for figure) 
fate_comp_combSubsets <- function(state1, state2) {
  
  # block on CD16 status
  print("combined CD16+/- analysis")
  subsetVar <- (dwnsmp$cluster == state1 | dwnsmp$cluster == state2) 
  modDS <- dwnsmp[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  comb <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"), block=paste(modDS$TP, modDS$sxCD16, sep="_"))
  comb[[1]]$gene <- rownames(comb[[1]])
  logFC <- findMarkers(dwn[,modDS$Barcode], modDS$cluster, test.type=c("t"), block=paste(modDS$TP, modDS$sxCD16, sep="_"))
  comb[[1]]$logFC <- logFC[[1]][match(comb[[1]]$gene, rownames(logFC[[1]])), 4]
  comb <- data.frame(comb[[1]][,c(5, 6, 4, 2, 3)])
  names(comb) <- c("gene",
                      paste(state1, state2, "logFC", sep="_"),
                      paste(state1, state2, "AUC", sep="_"),
                      paste(state1, state2, "pval", sep="_"), 
                      paste(state1, state2, "fdr", sep="_"))
  head(comb)
  
  return(comb)
}

comb.un.by <- fate_comp_combSubsets("unexposed", "bystander")
comb.un.inf <- fate_comp_combSubsets("unexposed", "infected")

test <- left_join(comb.un.by, comb.un.inf, by=c("gene"))
#*********************************************************************************** END (for figure)

# modify function for T4 high versus low IAV-transcribers
fate_comp_T4 <- function(state1, state2) {
  
  #CD16neg
  print("CD16neg subset")
  subsetVar <- (t4Inf$cluster == state1 | t4Inf$cluster == state2) & t4Inf$sxCD16=="FALSE"
  modDS <- t4Inf[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  CD16neg <- findMarkers(t4[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"))
  CD16neg[[1]]$gene <- rownames(CD16neg[[1]])
  logFC <- findMarkers(t4[,modDS$Barcode], modDS$cluster, test.type=c("t"))
  CD16neg[[1]]$logFC <- logFC[[1]][match(CD16neg[[1]]$gene, rownames(logFC[[1]])), 4]
  CD16neg <- data.frame(CD16neg[[1]][,c(5, 6, 4, 2, 3)])
  names(CD16neg) <- c("gene",
                      paste("CD16neg", state1, state2, "logFC", sep="_"),
                      paste("CD16neg", state1, state2, "AUC", sep="_"),
                      paste("CD16neg", state1, state2, "pval", sep="_"), 
                      paste("CD16neg", state1, state2, "fdr", sep="_"))
  head(CD16neg)
  
  #CD16pos
  print("CD16pos subset")
  subsetVar <- (t4Inf$cluster == state1 | t4Inf$cluster == state2) & t4Inf$sxCD16=="TRUE"
  modDS <- t4Inf[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  CD16pos <- findMarkers(t4[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"))
  CD16pos[[1]]$gene <- rownames(CD16pos[[1]])
  logFC <- findMarkers(t4[,modDS$Barcode], modDS$cluster, test.type=c("t"))
  CD16pos[[1]]$logFC <- logFC[[1]][match(CD16pos[[1]]$gene, rownames(logFC[[1]])), 4]
  CD16pos <- data.frame(CD16pos[[1]][,c(5, 6, 4, 2, 3)])
  names(CD16pos) <- c("gene",
                      paste("CD16pos", state1, state2, "logFC", sep="_"),
                      paste("CD16pos", state1, state2, "AUC", sep="_"),
                      paste("CD16pos", state1, state2, "pval", sep="_"), 
                      paste("CD16pos", state1, state2, "fdr", sep="_"))
  head(CD16pos)
  
  
  #Subset differences
  print("Testing for differences between subsets")
  subsetVar <- (t4Inf$cluster == state1 | t4Inf$cluster == state2)
  modDS <- t4Inf[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  dataList=list(apply(t4[,subsetVar], 1, function(y){
    mod <- lm(y ~ modDS$sxCD16 + modDS$cluster + modDS$sxCD16:modDS$cluster)
    c(summary(mod)$coefficients[4,4])}))
  df <- data.frame(dataList)
  df$gene <- rownames(df)
  df <- df[,c(2,1)]
  names(df) <- c("gene", "pval")
  df$fdr <- p.adjust(df$pval, method='fdr')
  names(df) <- c("gene", paste("subsetDiff", state1, state2, "pval", sep="_"), paste("subsetDiff", state1, state2, "fdr", sep="_"))
  head(df)
  
  ofInt <- full_join(CD16pos, CD16neg, by="gene")
  ofInt <- full_join(ofInt, df, by="gene")
  
  return(ofInt)
}

low.hi <- fate_comp_T4("low", "high")

# modify T4 high versus low IAV-transcribers function, combining subsets **************** (for figure)
fate_comp_T4_combSubsets <- function(state1, state2) {
  
  #block on CD16+/- statuse
  print("Time point 4")
  subsetVar <- (t4df$cluster == state1 | t4df$cluster == state2) #& t4Inf$sxCD16=="FALSE"
  modDS <- t4df[subsetVar,]
  modDS$cluster <- factor(modDS$cluster, levels=c(state1, state2))
  print(table(modDS$cluster))
  print(table(modDS$sxCD16))
  comb <- findMarkers(t4v2[,modDS$Barcode], modDS$cluster, test.type=c("wilcox"), block=modDS$sxCD16)
  comb[[1]]$gene <- rownames(comb[[1]])
  logFC <- findMarkers(t4v2[,modDS$Barcode], modDS$cluster, test.type=c("t"), block=modDS$sxCD16)
  comb[[1]]$logFC <- logFC[[1]][match(comb[[1]]$gene, rownames(logFC[[1]])), 4]
  comb <- data.frame(comb[[1]][,c(5, 6, 4, 2, 3)])
  names(comb) <- c("gene",
                      paste(state1, state2, "logFC", sep="_"),
                      paste(state1, state2, "AUC", sep="_"),
                      paste(state1, state2, "pval", sep="_"), 
                      paste(state1, state2, "fdr", sep="_"))
  head(comb)
  
  return(comb)
}

comb.low.hi <- fate_comp_T4_combSubsets("low", "high")
comb.low.by <- fate_comp_T4_combSubsets("low", "bystander")

test <- left_join(test, comb.low.by, by=c("gene"))
test <- left_join(test, comb.low.hi, by=c("gene"))

write.table(test, file=paste0(PopDiff, "results/subsampled_gene_averages_and_stats_for_22593hostgenes_blocked_by_subset.txt"), quote=F, row.names=F, col.names=T, qmethod="double", eol="\n", sep="\t", na="NA")
#*********************************************************************************** END (for figure)

#Combine stats with the averages
master <- full_join(master, un.by, by="gene")
master <- full_join(master, un.inf, by="gene")
master <- full_join(master, by.inf, by="gene")
master <- full_join(master, low.hi, by="gene")

#Add ensembl ids for genes
master$ens <- GA[match(master$gene, GA$geneAlt), 'stable']

#saveRDS(master, file=paste0(PopDiff, "results/subsampled_gene_averages_and_stats_for_22593hostgenes.rds"))
#master=readRDS(file=paste0(PopDiff, "results/subsampled_gene_averages_and_stats_for_22593hostgenes.rds"))

#Focus on highly expressed genes
master$expCut <- do.call(pmax, c(master[4:31], list(na.rm=TRUE))) 
subMaster <- dplyr::filter(master, expCut > 0.1)

####Unexposed vs bystander
#Recompute FDR & establish gene lists for unexposed vs bystander comparison
subMaster$CD16pos_unexposed_bystander_fdr <- p.adjust(subMaster$CD16pos_unexposed_bystander_pval, method="fdr")
subMaster$CD16neg_unexposed_bystander_fdr <- p.adjust(subMaster$CD16neg_unexposed_bystander_pval, method="fdr")
subMaster$subsetDiff_unexposed_bystander_fdr <- p.adjust(subMaster$subsetDiff_unexposed_bystander_pval, method="fdr")

# Add the fill variable: bystander, unexposed
subMaster$unexposed_bystander_fill <- ifelse((subMaster$CD16pos_unexposed_bystander_fdr >= 0.01 | subMaster$CD16neg_unexposed_bystander_fdr >= 0.01 | abs(subMaster$CD16pos_unexposed_bystander_logFC) < 0.2 | abs(subMaster$CD16neg_unexposed_bystander_logFC) < 0.2 | sign(subMaster$CD16pos_unexposed_bystander_logFC) != sign(subMaster$CD16neg_unexposed_bystander_logFC)), 'nonsig', ifelse(subMaster$CD16pos_unexposed_bystander_AUC < 0.5, 'bystander', 'unexposed'))

table(subMaster$unexposed_bystander_fill)

# Add the shape variable: CD16neg, CD16pos
subMaster$unexposed_bystander_shape <- ifelse(subMaster$subsetDiff_unexposed_bystander_fdr >= 0.01 | subMaster$unexposed_bystander_fill == "nonsig", "nonsig", ifelse(abs(subMaster$CD16neg_unexposed_bystander_logFC) > abs(subMaster$CD16pos_unexposed_bystander_logFC), 'CD16neg', 'CD16pos'))

table(subMaster$unexposed_bystander_shape)

#Define gene lists for GO enrichments
background=unique(subMaster$ens)

table(paste(subMaster$unexposed_bystander_fill, subMaster$unexposed_bystander_shape, sep="-"))

bystander = subMaster[subMaster$unexposed_bystander_fill=="bystander", 'ens']
unexposed = subMaster[subMaster$unexposed_bystander_fill=="unexposed", 'ens']
bystanderCD16neg = subMaster[subMaster$unexposed_bystander_fill=="bystander" & subMaster$unexposed_bystander_shape=="CD16neg", 'ens']
bystanderCD16pos = subMaster[subMaster$unexposed_bystander_fill=="bystander" & subMaster$unexposed_bystander_shape=="CD16pos", 'ens']
unexposedCD16neg = subMaster[subMaster$unexposed_bystander_fill=="unexposed" & subMaster$unexposed_bystander_shape=="CD16neg", 'ens']
unexposedCD16pos = subMaster[subMaster$unexposed_bystander_fill=="unexposed" & subMaster$unexposed_bystander_shape=="CD16pos", 'ens']

unexposed_bystander_enrLists = list(bystander, unexposed, bystanderCD16neg, bystanderCD16pos, unexposedCD16neg, unexposedCD16pos)
names(unexposed_bystander_enrLists) <- c("bystander", "unexposed", "bystander - CD16neg", "bystander - CD16pos", "unexposed - CD16neg", "unexposed - CD16pos")
str(unexposed_bystander_enrLists)

#Define function to perform GO enrichment
perform_GO_GR38 <- function(enrLists) {
  source(paste0(PopDiff, "/scripts/GOSeq_hg38.R"))
  library(data.table)
  paste("Performing enrichment")
  subGO=list()
  for (y in 1:length(enrLists)){
    print(y)
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
  print(str(subGO))
  res <- rbindlist(subGO, idcol=T, fill=T)
  return(res)
}

unexposed.bystander.GO.enrichment <- perform_GO_GR38(unexposed_bystander_enrLists)
table(unexposed.bystander.GO.enrichment$.id)

write.table(unexposed.bystander.GO.enrichment[,-10], paste0(PopDiff, "results/unexposed.bystander.logFC02.FDR01.GO.enrichment.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

####Unexposed vs infected
#Recompute FDR & establish gene lists for unexposed vs bystander comparison
subMaster$CD16pos_unexposed_infected_fdr <- p.adjust(subMaster$CD16pos_unexposed_infected_pval, method="fdr")
subMaster$CD16neg_unexposed_infected_fdr <- p.adjust(subMaster$CD16neg_unexposed_infected_pval, method="fdr")
subMaster$subsetDiff_unexposed_infected_fdr <- p.adjust(subMaster$subsetDiff_unexposed_infected_pval, method="fdr")

# Add the fill variable: infected, unexposed
subMaster$unexposed_infected_fill <- ifelse((subMaster$CD16pos_unexposed_infected_fdr >= 0.01 | subMaster$CD16neg_unexposed_infected_fdr >= 0.01 | abs(subMaster$CD16pos_unexposed_infected_logFC) < 0.2 | abs(subMaster$CD16neg_unexposed_infected_logFC) < 0.2 | sign(subMaster$CD16pos_unexposed_infected_logFC) != sign(subMaster$CD16neg_unexposed_infected_logFC)), 'nonsig', ifelse(subMaster$CD16pos_unexposed_infected_AUC < 0.5, 'infected', 'unexposed'))

table(subMaster$unexposed_infected_fill)

# Add the shape variable: CD16neg, CD16pos
subMaster$unexposed_infected_shape <- ifelse(subMaster$subsetDiff_unexposed_infected_fdr >= 0.01 | subMaster$unexposed_infected_fill == "nonsig", "nonsig", ifelse(abs(subMaster$CD16neg_unexposed_infected_logFC) > abs(subMaster$CD16pos_unexposed_infected_logFC), 'CD16neg', 'CD16pos'))

table(subMaster$unexposed_infected_shape)

#Define gene lists for GO enrichments
background=unique(subMaster$ens)

table(paste(subMaster$unexposed_infected_fill, subMaster$unexposed_infected_shape, sep="-"))

infected = subMaster[subMaster$unexposed_infected_fill=="infected", 'ens']
unexposed = subMaster[subMaster$unexposed_infected_fill=="unexposed", 'ens']
infectedCD16neg = subMaster[subMaster$unexposed_infected_fill=="infected" & subMaster$unexposed_infected_shape=="CD16neg", 'ens']
infectedCD16pos = subMaster[subMaster$unexposed_infected_fill=="infected" & subMaster$unexposed_infected_shape=="CD16pos", 'ens']
unexposedCD16neg = subMaster[subMaster$unexposed_infected_fill=="unexposed" & subMaster$unexposed_infected_shape=="CD16neg", 'ens']
unexposedCD16pos = subMaster[subMaster$unexposed_infected_fill=="unexposed" & subMaster$unexposed_infected_shape=="CD16pos", 'ens']


unexposed_infected_enrLists = list(infected, unexposed, infectedCD16neg, infectedCD16pos, unexposedCD16neg, unexposedCD16pos)
names(unexposed_infected_enrLists) <- c("infected", "unexposed", "infected - CD16neg", "infected - CD16pos", "unexposed - CD16neg", "unexposed - CD16pos")
str(unexposed_infected_enrLists)

unexposed.infected.GO.enrichment <- perform_GO_GR38(unexposed_infected_enrLists)
table(unexposed.infected.GO.enrichment$.id)

write.table(unexposed.infected.GO.enrichment[,-10], paste0(PopDiff, "results/unexposed.infected.logFC02.FDR01.GO.enrichment.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

####bystander vs infected
#Recompute FDR & establish gene lists for bystander vs bystander comparison
subMaster$CD16pos_bystander_infected_fdr <- p.adjust(subMaster$CD16pos_bystander_infected_pval, method="fdr")
subMaster$CD16neg_bystander_infected_fdr <- p.adjust(subMaster$CD16neg_bystander_infected_pval, method="fdr")
subMaster$subsetDiff_bystander_infected_fdr <- p.adjust(subMaster$subsetDiff_bystander_infected_pval, method="fdr")

# Add the fill variable: infected, bystander
subMaster$bystander_infected_fill <- ifelse((subMaster$CD16pos_bystander_infected_fdr >= 0.01 | subMaster$CD16neg_bystander_infected_fdr >= 0.01 | abs(subMaster$CD16pos_bystander_infected_logFC) < 0.2 | abs(subMaster$CD16neg_bystander_infected_logFC) < 0.2 | sign(subMaster$CD16pos_bystander_infected_logFC) != sign(subMaster$CD16neg_bystander_infected_logFC)), 'nonsig', ifelse(subMaster$CD16pos_bystander_infected_AUC < 0.5, 'infected', 'bystander'))

table(subMaster$bystander_infected_fill)

# Add the shape variable: CD16neg, CD16pos
subMaster$bystander_infected_shape <- ifelse(subMaster$subsetDiff_bystander_infected_fdr >= 0.01 | subMaster$bystander_infected_fill == "nonsig", "nonsig", ifelse(abs(subMaster$CD16neg_bystander_infected_logFC) > abs(subMaster$CD16pos_bystander_infected_logFC), 'CD16neg', 'CD16pos'))

table(subMaster$bystander_infected_shape)

#Define gene lists for GO enrichments
background=unique(subMaster$ens)

table(paste(subMaster$bystander_infected_fill, subMaster$bystander_infected_shape, sep="-"))

infected = subMaster[subMaster$bystander_infected_fill=="infected", 'ens']
bystander = subMaster[subMaster$bystander_infected_fill=="bystander", 'ens']
infectedCD16neg = subMaster[subMaster$bystander_infected_fill=="infected" & subMaster$bystander_infected_shape=="CD16neg", 'ens']
infectedCD16pos = subMaster[subMaster$bystander_infected_fill=="infected" & subMaster$bystander_infected_shape=="CD16pos", 'ens']
bystanderCD16neg = subMaster[subMaster$bystander_infected_fill=="bystander" & subMaster$bystander_infected_shape=="CD16neg", 'ens']
bystanderCD16pos = subMaster[subMaster$bystander_infected_fill=="bystander" & subMaster$bystander_infected_shape=="CD16pos", 'ens']


bystander_infected_enrLists = list(infected, bystander, infectedCD16neg, infectedCD16pos, bystanderCD16neg, bystanderCD16pos)
names(bystander_infected_enrLists) <- c("infected", "bystander", "infected - CD16neg", "infected - CD16pos", "bystander - CD16neg", "bystander - CD16pos")
str(bystander_infected_enrLists)

bystander.infected.GO.enrichment <- perform_GO_GR38(bystander_infected_enrLists)
table(bystander.infected.GO.enrichment$.id)

write.table(bystander.infected.GO.enrichment[,-10], paste0(PopDiff, "results/bystander.infected.logFC02.FDR01.GO.enrichment.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

####low vs high
#Recompute FDR & establish gene lists for low vs low comparison
subMaster$CD16pos_low_high_fdr <- p.adjust(subMaster$CD16pos_low_high_pval, method="fdr")
subMaster$CD16neg_low_high_fdr <- p.adjust(subMaster$CD16neg_low_high_pval, method="fdr")
subMaster$subsetDiff_low_high_fdr <- p.adjust(subMaster$subsetDiff_low_high_pval, method="fdr")

# Add the fill variable: high, low
subMaster$low_high_fill <- ifelse((subMaster$CD16pos_low_high_fdr >= 0.01 | subMaster$CD16neg_low_high_fdr >= 0.01 | abs(subMaster$CD16pos_low_high_logFC) < 0.2 | abs(subMaster$CD16neg_low_high_logFC) < 0.2 | sign(subMaster$CD16pos_low_high_logFC) != sign(subMaster$CD16neg_low_high_logFC)), 'nonsig', ifelse(subMaster$CD16pos_low_high_AUC < 0.5, 'high', 'low'))

table(subMaster$low_high_fill)

# Add the shape variable: CD16neg, CD16pos
subMaster$low_high_shape <- ifelse(subMaster$subsetDiff_low_high_fdr >= 0.01 | subMaster$low_high_fill == "nonsig", "nonsig", ifelse(abs(subMaster$CD16neg_low_high_logFC) > abs(subMaster$CD16pos_low_high_logFC), 'CD16neg', 'CD16pos'))

table(subMaster$low_high_shape)

#Define gene lists for GO enrichments
background=unique(subMaster$ens)

table(paste(subMaster$low_high_fill, subMaster$low_high_shape, sep="-"))

high = subMaster[subMaster$low_high_fill=="high", 'ens']
low = subMaster[subMaster$low_high_fill=="low", 'ens']
highCD16neg = subMaster[subMaster$low_high_fill=="high" & subMaster$low_high_shape=="CD16neg", 'ens']
highCD16pos = subMaster[subMaster$low_high_fill=="high" & subMaster$low_high_shape=="CD16pos", 'ens']
lowCD16neg = subMaster[subMaster$low_high_fill=="low" & subMaster$low_high_shape=="CD16neg", 'ens']
lowCD16pos = subMaster[subMaster$low_high_fill=="low" & subMaster$low_high_shape=="CD16pos", 'ens']


low_high_enrLists = list(high, low, highCD16neg, highCD16pos, lowCD16neg, lowCD16pos)
names(low_high_enrLists) <- c("high", "low", "high - CD16neg", "high - CD16pos", "low - CD16neg", "low - CD16pos")
str(low_high_enrLists)

low.high.GO.enrichment <- perform_GO_GR38(low_high_enrLists)
table(low.high.GO.enrichment$.id)

write.table(low.high.GO.enrichment[,-10], paste0(PopDiff, "results/low.high.logFC02.FDR01.GO.enrichment.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#save the final dataframe
write.table(subMaster, file=paste0(PopDiff, "results/subsampled_gene_averages_and_stats_for_6669hostgenes.txt"), quote=F, row.names=F, col.names=T, qmethod="double", eol="\n", sep="\t", na="NA")

#write gene lists to file
write.table(background, file=paste0(PopDiff, "results/gene_lists/subsampled/background_6669genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")

write.table(bystander_infected_enrLists[["infected"]], file=paste0(PopDiff, "results/gene_lists/subsampled/infected_relative_to_bystander_155genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")
write.table(bystander_infected_enrLists[["bystander"]], file=paste0(PopDiff, "results/gene_lists/subsampled/bystander_relative_to_infected_137genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")

write.table(unexposed_bystander_enrLists[["unexposed"]], file=paste0(PopDiff, "results/gene_lists/subsampled/unexposed_relative_to_bystander_118genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")
write.table(unexposed_bystander_enrLists[["bystander"]], file=paste0(PopDiff, "results/gene_lists/subsampled/bystander_relative_to_unexposed_523genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")

write.table(unexposed_infected_enrLists[["unexposed"]], file=paste0(PopDiff, "results/gene_lists/subsampled/unexposed_relative_to_infected_144genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")
write.table(unexposed_infected_enrLists[["infected"]], file=paste0(PopDiff, "results/gene_lists/subsampled/infected_relative_to_unexposed_515genes.txt"), quote=F, row.names=F, col.names=F, eol="\n", sep="\n", na="NA")

##################################################################################################################
############ Focus on the subset-specific responses to IAV  ######################################################
##################################################################################################################
#Filter for those genes that exhibit a subset-specific response
mOpp <-filter(subMaster, subsetDiff_unexposed_bystander_fdr < 0.01 | subsetDiff_unexposed_infected_fdr < 0.01) 
rownames(mOpp) <- mOpp$gene
mOpp <- mOpp[,2:27]
head(mOpp)

#Isolate the 335 genes 
opp <- dwn[rownames(mOpp),]
dim(opp)

#perform clustering on scaled expression values
dSub <- dist(t(scale(t(opp))))
hcSub <- hclust(dSub, method="ward.D2")
plot(hcSub)

clustSub <- cutreeDynamic(hcSub, distM=as.matrix(dSub),
                          minClusterSize=2, deepSplit=1)
table(clustSub)

df <- mOpp
df$cluster <- clustSub

#Plot the average gene expression patterns
toLab <- c("CH25H","CCL7","CCL8","CXCL10","LYZ","CCL13","ISG20","TNF","IL6","STAT1","ITIF1B","IFIT5","IFI35","ISG15","ISG20","APOBEC3G","S100A9","S100A12","S100A8","CX3CR1", "CEBPA", "JUN", "EGR2", "SIRT1", "MS4A4A")
labs = data.frame(gene=as.character(rownames(mOpp)), idx=1:335)
tolabDF = filter(labs, gene %in% toLab)

top_genes = rowAnnotation(idx = anno_mark(at = tolabDF$idx, labels = tolabDF$gene))

hm1 <- Heatmap(t(scale(t(as.matrix(mOpp)))), 
               cluster_rows=hcSub,
               cluster_columns = F,
               column_split = factor(paste(sapply(strsplit(sapply(strsplit(as.character(colnames(mOpp)), "_"), "[[", 3), "_"), "[[", 1),
                                           sapply(strsplit(sapply(strsplit(as.character(colnames(mOpp)), "_"), "[[", 2), "_"), "[[", 1),
                                           sep="_"), levels=c("CD16neg_T0", "CD16neg_unexposed", "CD16neg_bystander", "CD16neg_infected",
                                                              "CD16pos_T0", "CD16pos_unexposed", "CD16pos_bystander", "CD16pos_infected")),
               right_annotation = rowAnnotation(cluster=as.factor(clustSub), idx = anno_mark(at = tolabDF$idx, labels = tolabDF$gene)),
               show_column_names = F,
               show_row_names=F,
               col = inlmisc::GetColors(256, scheme = "BuRd"),
               heatmap_legend_param = list(
               title = "scaled", at = c(-6, 0, 6))
)
hm1

#Extract info for supplementary table
tabS3C <- mOpp
tabS3C$module_pre_reordering <- clustSub
tabS3C$bystander <- ifelse(subMaster[match(rownames(mOpp), subMaster$gene), 'subsetDiff_unexposed_bystander_fdr'] < 0.01, TRUE, FALSE)
tabS3C$infected <- ifelse(subMaster[match(rownames(mOpp), subMaster$gene), 'subsetDiff_unexposed_infected_fdr'] < 0.01, TRUE, FALSE)
table(paste(tabS3C$infected, tabS3C$bystander))
tabS3C <- tabS3C[hcSub$labels[hcSub$order],]

write.table(tabS3C, file=paste0(PopDiff, "results/subset-specific_genes_335genes.txt"), row.names=T, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#Perform GO enrichments on the modules
background=unique(subMaster$ens)

tabS3C$ens = GA[match(rownames(tabS3C), GA$geneAlt), 'stable']

clustList=list()
for (i in 1:8) {
  clustList[[i]] = tabS3C[tabS3C$module_pre_reordering==i, 'ens']
}

#Define function to perform GO enrichment
perform_GO_GR38 <- function(enrLists) {
  source(paste0(PopDiff, "scripts/GOSeq_hg38.R"))
  library(data.table)
  paste("Performing enrichment")
  subGO=list()
  for (y in 1:length(enrLists)){
    print(y)
    t=unique(enrLists[[y]])
    goRES <- GOSeq(t, background, addCI=TRUE)
    print(str(goRES))
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
    nm=y
    subGO[[nm]]=goRES
    
  }
  print(str(subGO))
  res <- rbindlist(subGO, idcol=T, fill=T)
  return(res)
}

subsetDiff.enrichment <- perform_GO_GR38(clustList)
table(subsetDiff.enrichment$.id)

write.table(subsetDiff.enrichment[,-10], paste0(PopDiff, "results/subsetDiff.GO.enrichment.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

##################################################################################################################
############ End  ################################################################################################
##################################################################################################################