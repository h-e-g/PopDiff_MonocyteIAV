
get_proportions = function(rawCounts_DT, bulk_TPM_matrix, cell_annotation, sample_ids, normalize=TRUE, log_transform = TRUE){
  # rawCounts_DT: data table containing the raw count for each gene x barcode
  # bulk_TPM_matrix: matrix of Bulk TPM from samples on which we want to estimate the proportions
  #cell_annotation : Detail of assignment of each cell to a cluster/individual/condition/TP
  #sample_ids vector of the form individuals to consider when defining pseudo-Bulk profiles in each cluster

  get_clusterTPM=function(rawCounts_DT, sample_ids){
    # compute pseudoBulk profile for each cluster
    clust_TPM=rawCounts_DT[ID%in%sample_ids & cluster!='Lymphocytes',]
    setkey(clust_TPM,gene, cluster)
    clust_TPM=clust_TPM[,.(count=sum(count)),by=.(gene, cluster)]
    clust_TPM[,TPM:=count/sum(count)*1e6,by=.(cluster)]
    # clust_TPM[,count_allCluster:=sum(count)*1e6,by=.(TP,cluster)]
    clust_TPM
  }
  clust_TPM=get_clusterTPM(rawCounts_DT, sample_ids)

  get_clusterWeightsReads = function(sample_ids, cell_annotation){
      # compute total number of UMIs in a cluster, for selected individual
      counts=cell_annotation[,.(nbReads = sum(total_counts)),by = .(ID,cluster)]
      # keep only selected samples, and compute total number of UMIs in that cluster
      counts = counts[ID %in% sample_ids, .(nbReads = sum(nbReads)), by=cluster]
      # trick to force each cluster to be present in the final table, even if there are no cells for this cluster
      cluster_levels=unique(counts[,cluster])
      counts=merge(counts, data.table(cluster = cluster_levels, nbReads = 0) ,all=T, by=c('cluster','nbReads')) [, .(count = sum(nbReads) ), by=cluster]
      # compute cluster weight as % of all transcripts that are in that cluster
      counts[,weight:= count/sum(count)]
      counts[is.na(cluster),cluster:='NA']
      counts
  }

  clusterWeights=get_clusterWeightsReads(sample_ids, cell_annotation)

	  #clusterLevels=clust_TPM[,unique(cluster)]
    # compute global pseudoBulk estimate from estiamtes of each cluster
    clust_TPM=merge(clust_TPM,clusterWeights[,.(cluster,weight)],by='cluster')
    clust_TPM[,TPM_pseudoBulk:=sum(TPM*weight),by=gene]

    # get pseudoBulk estimates per cluster and overall
    Cell_TPM = dcast(clust_TPM, gene+TPM_pseudoBulk~cluster,value.var='TPM',fill=0)

    # make sure that genes are in the same order
    study_genes=intersect(Cell_TPM$gene,rownames(bulk_TPM_matrix))
    bulk_TPM_matrix=bulk_TPM_matrix[match(study_genes,rownames(bulk_TPM_matrix)),]
    Cell_TPM = Cell_TPM[match(study_genes,Cell_TPM$gene),]

    clusterLevels=clust_TPM[,unique(cluster)]

    # extract cluster TPMs
    Clust_TPM_mat = as.data.frame( Cell_TPM[, mget(clusterLevels) ] )
    Clust_TPM_mat = as.matrix(Clust_TPM_mat)

    log2_transform = function(x){ log2(1+x) }
    log2_untransform = function(x){ pmax(2^x-1 , 0) }

    TPM_pseudoBulk = Cell_TPM$TPM_pseudoBulk
    if(log_transform){
        TPM_pseudoBulk = log2_transform( TPM_pseudoBulk )
        bulk_TPM_matrix = log2_transform( bulk_TPM_matrix )
        Clust_TPM_mat = log2_transform( Clust_TPM_mat )
    }
    mean_TPM_bulk = apply(bulk_TPM_matrix,1,mean)
    if( normalize ){
        gene_scaling = mean_TPM_bulk - TPM_pseudoBulk
        TPM_bulk_scaled = bulk_TPM_matrix - gene_scaling
    }else{
        TPM_bulk_scaled = bulk_TPM_matrix
    }

    # run DeconRNASeq
    DF = as.data.frame(TPM_bulk_scaled)
    rownames(DF) = study_genes
	  Clust_TPM_mat = as.data.frame(Clust_TPM_mat)
    rownames(Clust_TPM_mat) = study_genes
    result = DeconRNASeq(DF,Clust_TPM_mat)
    proportions = result$out.all
    rownames(proportions) = colnames(bulk_TPM_matrix)
    colnames(proportions) = colnames(Clust_TPM_mat)
    return( proportions )
}

CellMetaData=fread('./data/CellMetaData.txt',sep='\t')
rawCounts_Cells=fread('gunzip -c ./data/rawCounts_Cells.txt.gz',sep='\t')
bulk_TPM=fread('gunzip -c ./data/bulk_TPM.txt.gz',sep='\t')
bulk_TPM_mat=as.matrix(bulk_TPM[,-1])
rownames(bulk_TPM_mat)=bulk_TPM$ID

high_responders=c("EUB058", "EUB078", "AFB186", "AFB069")
low_responders=c("EUB094", "AFB108", "EUB061", "AFB105")
sample_ids=c(low_responders,high_responders)

library(DeconRNASeq)

proportions_IAV = get_proportions(rawCounts_Cells, bulk_TPM_mat, CellMetaData, sample_ids, normalize=TRUE, log_transform = TRUE)

fwrite(data.table(ID=rownames(proportions_IAV),proportions_IAV),file='./results/CellStateProportions.txt',sep='\t')
