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

#Load miscellaneous other relevant information
GA=readRDS(paste0(PopDiff, "data/gene_annotation_22593hostgenes.rds"))
flu <- read.table(paste0(PopDiff, "data/EvoImmunoPop_perFlu_estimates.txt"), header=T)
t0.master=read.table(paste0(PopDiff, "results/T0_three_subsets_stats_4589genes.compare.txt"), header=T, na.strings="NA", stringsAsFactors = F) #generated with script 01

##################################################################################################################
############ T0 inter-individual variability  ####################################################################
##################################################################################################################
#Extract T0
t0 <- filter(keepMeta, TP==0)
t0mat <- mat[t0.master$gene,t0$Barcode]

################
################ T0 individual differences in proportions of subsets
################
#calculate the number of cells per individual at T0
cnts <- data.frame(table(t0$ID1))

#calculate the number and proportion of each of the three subsets for each individual at T0
comp <- t0 %>% group_by(ID1) %>% dplyr::count(T0subset) %>% data.frame()
comp$total <- cnts[match(comp$ID1, cnts$Var1),'Freq']
comp$per <- comp$n/comp$total

#look for correlation to flu phenotypes
compDF <- dcast(comp[,c(1, 2, 5)], ID1 ~ T0subset)
compDF$EIP <- flu[match(paste(compDF$ID1, "5", sep="-"), flu$sample), 'perUniq']

cor.test(compDF$classical, compDF$EIP)
cor.test(compDF$nonclassical, compDF$EIP)
cor.test(compDF$intermediate, compDF$EIP) #intermediate monocytes again correlate with flu levels

perInf <- filter(keepMeta, TP != 0 & COND=="IAV") %>% group_by(ID1, TP) %>% dplyr::count(sxFlu) %>% data.frame() %>% dcast(ID1+TP~sxFlu)
perInf$inf <- perInf$`TRUE`/(perInf$`TRUE`+perInf$`FALSE`)
perInf <- dcast(perInf[,c(1,2,5)], ID1~TP)

identical(compDF$ID1, perInf$ID1)
compDF <- cbind(compDF, perInf)

#look per time-point
cor.test(compDF$intermediate, compDF$`2`)
cor.test(compDF$intermediate, compDF$`4`) #**only significant 
cor.test(compDF$intermediate, compDF$`6`)
cor.test(compDF$intermediate, compDF$`8`)

cor.test(compDF$classical, compDF$`2`)
cor.test(compDF$classical, compDF$`4`) 
cor.test(compDF$classical, compDF$`6`)
cor.test(compDF$classical, compDF$`8`)

cor.test(compDF$nonclassical, compDF$`2`)
cor.test(compDF$nonclassical, compDF$`4`) 
cor.test(compDF$nonclassical, compDF$`6`)
cor.test(compDF$nonclassical, compDF$`8`)

perHigh <- filter(keepMeta, TP == 4 & COND=="IAV") %>% group_by(ID1) %>% dplyr::count(cluster) %>% data.frame() %>% dcast(ID1~cluster)
perHigh$inf <- (perHigh$low+ perHigh$high)/(perHigh$low+perHigh$high+perHigh$bystander)
perHigh$perHigh <- (perHigh$high)/(perHigh$low+perHigh$high)

cor.test(compDF$nonclassical, perHigh$perHigh)
cor.test(compDF$classical, perHigh$perHigh)
cor.test(compDF$intermediate, perHigh$perHigh)

################
################ T0 individual differences in gene expression
################
# Perform Kruskal-Wallis Rank Sum Test to identify genes differently expressed between individuals
classical=list(apply(t0mat[,t0$T0subset=="classical"], 1, function(y){
  kt <- kruskal.test(y, t0[t0$T0subset=="classical", 'ID1'])
  c(kt$p.value, kt$statistic)}))
classical <- data.frame(t(data.frame(classical)))
classical$fdr <- p.adjust(classical$V1, method="fdr")
table(classical$fdr<0.01)

intermediate=list(apply(t0mat[,t0$T0subset=="intermediate"], 1, function(y){
  kt <- kruskal.test(y, t0[t0$T0subset=="intermediate", 'ID1'])
  c(kt$p.value, kt$statistic)}))
intermediate <- data.frame(t(data.frame(intermediate)))
intermediate$fdr <- p.adjust(intermediate$V1, method="fdr")
table(intermediate$fdr<0.05)

nonclassical=list(apply(t0mat[,t0$T0subset=="nonclassical"], 1, function(y){
  kt <- kruskal.test(y, t0[t0$T0subset=="nonclassical", 'ID1'])
  c(kt$p.value, kt$statistic)}))
nonclassical <- data.frame(t(data.frame(nonclassical)))
nonclassical$fdr <- p.adjust(nonclassical$V1, method="fdr")
table(nonclassical$fdr<0.01)

#Find the averages
tmpList = list()
for (i in levels(t0$ID1)){
  avg <- unname(rowMeans(t0mat[, t0$T0subset=="classical" & t0$ID1==i]))
  tmpList[[i]] = avg
}
tmpDF <- data.frame(tmpList)
tmpDF$gene <- rownames(t0mat)
classical <- cbind(classical, tmpDF)

tmpList = list()
for (i in levels(t0$ID1)){
  avg <- unname(rowMeans(t0mat[, t0$T0subset=="intermediate" & t0$ID1==i]))
  tmpList[[i]] = avg
}
tmpDF <- data.frame(tmpList)
tmpDF$gene <- rownames(t0mat)
intermediate <- cbind(intermediate, tmpDF)

tmpList = list()
for (i in levels(t0$ID1)){
  avg <- unname(rowMeans(t0mat[, t0$T0subset=="nonclassical" & t0$ID1==i]))
  tmpList[[i]] = avg
}
tmpDF <- data.frame(tmpList)
tmpDF$gene <- rownames(t0mat)
nonclassical <- cbind(nonclassical, tmpDF)

#add correlation to flu
phenotype <- flu[match(paste0(levels(t0$ID1), "-5"), flu$sample),'perUniq']

classical$fluCor_pval=apply(classical[,4:11], 1, function(y){
  mod=summary(lm(y~phenotype));
  c(mod$coeff[2,4])})

classical$fluCor_coeff=apply(classical[,4:11], 1, function(y){
  mod=cor.test(y, phenotype);
  c(mod$estimate)})

intermediate$fluCor_pval=apply(intermediate[,4:11], 1, function(y){
  mod=summary(lm(y~phenotype));
  c(mod$coeff[2,4])})

intermediate$fluCor_coeff=apply(intermediate[,4:11], 1, function(y){
  mod=cor.test(y, phenotype);
  c(mod$estimate)})

nonclassical$fluCor_pval=apply(nonclassical[,4:11], 1, function(y){
  mod=summary(lm(y~phenotype));
  c(mod$coeff[2,4])})

nonclassical$fluCor_coeff=apply(nonclassical[,4:11], 1, function(y){
  mod=cor.test(y, phenotype);
  c(mod$estimate)})

#Add the ensemble ids
classical$ens <- GA[match(classical$gene, GA$geneAlt), 'stable']
intermediate$ens <- GA[match(intermediate$gene, GA$geneAlt), 'stable']
nonclassical$ens <- GA[match(nonclassical$gene, GA$geneAlt), 'stable']

#filter
classical2 <- filter(classical, fdr < 0.01)
intermediate2 <- filter(intermediate, fdr < 0.01)
nonclassical2 <- filter(nonclassical, fdr < 0.01)

classical2$fluCor_fdr <- p.adjust(classical2$fluCor_pval, method="fdr")
intermediate2$fluCor_fdr <- p.adjust(intermediate2$fluCor_pval, method="fdr")
nonclassical2$fluCor_fdr <- p.adjust(nonclassical2$fluCor_pval, method="fdr")

arrange(nonclassical2, fluCor_pval) %>% head()
arrange(intermediate2, fluCor_pval) %>% head()
arrange(classical2, fluCor_pval) %>% head()

table(classical2$fluCor_pval < 0.01) #126
table(intermediate2$fluCor_pval < 0.01) #6
table(nonclassical2$fluCor_pval < 0.01) #18

#Add more information to the the master T0 table
t0.master$classical_id_pval <- classical[match(t0.master$ens, classical$ens), 'V1']
t0.master$classical_id_fdr <- classical[match(t0.master$ens, classical$ens), 'fdr']
t0.master$classical_fluCor_pval <- classical[match(t0.master$ens, classical$ens), 'fluCor_pval']
t0.master$classical_fluCor_cc <- classical[match(t0.master$ens, classical$ens), 'fluCor_coeff']

t0.master$intermediate_id_pval <- intermediate[match(t0.master$ens, intermediate$ens), 'V1']
t0.master$intermediate_id_fdr <- intermediate[match(t0.master$ens, intermediate$ens), 'fdr']
t0.master$intermediate_fluCor_pval <- intermediate[match(t0.master$ens, intermediate$ens), 'fluCor_pval']
t0.master$intermediate_fluCor_cc <- intermediate[match(t0.master$ens, intermediate$ens), 'fluCor_coeff']

t0.master$nonclassical_id_pval <- nonclassical[match(t0.master$ens, nonclassical$ens), 'V1']
t0.master$nonclassical_id_fdr <- nonclassical[match(t0.master$ens, nonclassical$ens), 'fdr']
t0.master$nonclassical_fluCor_pval <- nonclassical[match(t0.master$ens, nonclassical$ens), 'fluCor_pval']
t0.master$nonclassical_fluCor_cc <- nonclassical[match(t0.master$ens, nonclassical$ens), 'fluCor_coeff']

classical$mean <- rowMeans(classical[,4:11], na.rm=T)
classical$sd <- rowSds(as.matrix(classical[,4:11]), na.rm=T)
intermediate$mean <- rowMeans(intermediate[,4:11], na.rm=T)
intermediate$sd <- rowSds(as.matrix(intermediate[,4:11]), na.rm=T)
nonclassical$mean <- rowMeans(nonclassical[,4:11], na.rm=T)
nonclassical$sd <- rowSds(as.matrix(nonclassical[,4:11]), na.rm=T)

t0.master$classical_mean <- classical[match(t0.master$ens, classical$ens), 'mean']
t0.master$classical_sd <- classical[match(t0.master$ens, classical$ens), 'sd']
t0.master$intermediate_mean <- intermediate[match(t0.master$ens, intermediate$ens), 'mean']
t0.master$intermediate_sd <- intermediate[match(t0.master$ens, intermediate$ens), 'sd']
t0.master$nonclassical_mean <- nonclassical[match(t0.master$ens, nonclassical$ens), 'mean']
t0.master$nonclassical_sd <- nonclassical[match(t0.master$ens, nonclassical$ens), 'sd']

write.table(t0.master, file=paste0(PopDiff, "results/T0_three_subsets_stats_4589genes_wInterindiviudalVariability.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

#Generate candidate gene lists
#Define lists for enrichments
nonclassicList <- nonclassical2[nonclassical2$fluCor_pval < 0.01, 'ens']
classicList <- classical2[classical2$fluCor_pval < 0.01, 'ens']
intermediateList <- intermediate2[intermediate2$fluCor_pval < 0.01, 'ens']

interest <- unique(c(enrLists[[1]], enrLists[[2]], enrLists[[3]])) #135 genes sig diff between donors & nominally correlates with flu in at least one subset

source(paste0(PopDiff, "/scripts/GOSeq_hg38.R"))

background = t0.master$ens

goRES <- GOSeq(interest, background, addCI=TRUE)
goRES$genes <- NA
tmpGA <- dplyr::filter(GA, stable %in% interest)
for (i in 1:length(goRES$category)) {
  category=goRES$category[i]
  genes <- tmpGA[grepl(category, tmpGA$GO),'Symbol']
  genes <- paste(genes, collapse=",")
  goRES[i, 'genes'] = genes
}

write.table(goRES[,-9], paste0(PopDiff, "results/T0_fluCor_subsetspecific_GO.enrichments_4589backgroundgenes_ONELIST.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))

interestGA <- filter(GA, stable %in% interest)
keep <- interestGA[grepl("GO:0051607|GO:0060337", interestGA$GO),'Symbol']

table1 <- filter(t0.master, gene %in% keep)
table1 <- table1[,c(2,1,7,16:17,19:20,22:23,25:32)]
table1 <- table1[,c(1:5,12,13,6,7,14,15,8,9,16,17,10,11)]
table1

table1$minP <- do.call(pmin, table1[c(9,13,17)])
table1 <- arrange(table1, minP)

write.table(table1, file=paste0(PopDiff, "results/table1_19genes_raw.txt"), row.names=F, col.names=T, sep="\t", eol="\n", quote=F, qmethod=c("double"))
write.table(interest, file=paste0(PopDiff, "results/135gene_candidates_ensemblID.txt"), row.names=F, col.names=F, sep="\n", eol="\n", quote=F, qmethod=c("double"))

#Plot IFITM3
t0$T0subset <- factor(t0$T0subset, levels=c("classical", "intermediate", "nonclassical"))
IFITM3 <- ggplot(t0, aes(x=T0subset, y=t0mat["IFITM3",])) +
  geom_violin(aes(fill=ID1, color=ID1), scale="width") +
  #geom_boxplot(aes(fill=ID1, color=ID1), notch=T) +
  theme_light() +
  theme(panel.grid=element_blank(),
        legend.position='none') +
  xlab("") +
  ylab("IFITM3") +
  scale_color_viridis(direction=-1, option="C", discrete=T, end=0.9) +
  scale_fill_viridis(direction=-1, option="C", discrete=T, end=0.9)
IFITM3

##################################################################################################################
############ End #################################################################################################
##################################################################################################################
