
library(data.table)

regulonTargetsInfo=fread('./data/regulonTargetsInfo.txt')

background=fread('./results/subsampled_gene_averages_and_stats_for_6669hostgenes.compare.txt')
background=background$gene
subsetSpecific_genes=fread('./results/subset-specific_genes_335genes.compare.txt')
colnames(subsetSpecific_genes)[1]='gene'

odds.ratio=function(tab, alpha = 0.05){
 test=fisher.test(tab,conf.level=1-alpha)
  oframe <- data.frame(LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], alpha = alpha,P=test$p.value)
  oframe
}

TestList=function(TFtargets,geneList,background){
  # geneList=intersect(geneList,regulonTargetsInfo[ highConfAnnot==TRUE , gene])
  # background=intersect(background,regulonTargetsInfo[ highConfAnnot==TRUE, gene])
  tab=table(background%in%geneList,background%in%TFtargets)
  data.table(PctBinding=mean(geneList%in%TFtargets),PctBindingBackground=mean(background%in%TFtargets),odds.ratio(tab))
    }

Enrich_clust=list()
  for (i in 1:8){
    geneList=subsetSpecific_genes[module_pre_reordering==i,gene]
    Enrich_clust[[i]]=regulonTargetsInfo[ highConfAnnot==TRUE ,TestList(gene,geneList,background),by=TF][order(-LowerCI)]
    Enrich_clust[[i]][,module_pre_reordering:=i]
    cat('cluster',i,'done\n')
    }
Enrich_clust=rbindlist(Enrich_clust)
Enrich_clust[,FDR:=p.adjust(P,'fdr')]
fwrite(Enrich_clust[FDR<.01,],file='./results/TF_enrichment.subset_differential_clusters.txt',sep='\t')

T0_background=fread('./results/T0_three_subsets_stats_4589genes_wInterindiviudalVariability.compare.txt')
T0_candidates=fread('./results/135gene_candidates_ensemblID.txt',header=F)$V1
T0_background_symbol=T0_background$gene
T0_candidates_symbol=T0_background$gene[match(T0_candidates,T0_background$ens)]

Enrich_T0=regulonTargetsInfo[ highConfAnnot==TRUE, TestList(gene, T0_candidates_symbol, T0_background_symbol), by=TF ][order( -LowerCI )]
Enrich_T0[, FDR:=p.adjust(P,'fdr')]
fwrite(Enrich_T0, file='./results/T0_TF.enrichment.txt',sep='\t')
