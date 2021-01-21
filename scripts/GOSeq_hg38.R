GOSeq=function(geneList,background,FDR=0.05,biasCorrect=TRUE,overOnly=T,addCI=F,bias.data=NULL,...){
  # geneList : vector with Ensembl gene ID of target genes
  # background : vector with Ensembl gene ID of background genes
  # FDR: threshold for adjusted enrichment Pvalue to be considered significant
  # biasCorrect: logical, if TRUE, enrichment will be adjusted on gene length or the variable provided in bias.data. default TRUE
  # overOnly: logical,  if TRUE, only over representation is tested. else underrepresentation is tested too. default TRUE
  #### addGenes: should gene participating to GO enrichment be reported. default TRUE. ****removed for compatibility
  # addCI: should odds ratio and associated confidence interval be reported. default FALSE .
  # bias.data : variable to adjust the enrichment on.
  require(GO.db)
  require(goseq)
  require(data.table)
  
  odds.ratio=function(tab, alpha = 0.05){
    test=fisher.test(tab,conf.level=1-alpha)
    oframe <- data.frame(LowerCI = test$conf.int[1], OR = test$est, UpperCI = test$conf.int[2], alpha = alpha,P=test$p.value)
    oframe
  }
  counts2tab=function(xy,x,y,tot){matrix(c(tot-x-y+xy,x-xy,y-xy,xy),2)}
  
  if(is.numeric(geneList) | is.logical(geneList)){
    geneList=background[geneList]
  }
  DE=background%in%geneList
  names(DE)=background
  if(is.null(bias.data)){
    nullP=nullp(DE, 'hg38', 'ensGene',plot.fit=TRUE)
  }else{
    nullP=nullp(DE, 'hg38', 'ensGene',plot.fit=TRUE,...)
  }
  if(!biasCorrect){
    nullP$pwf=mean(DE)
  }
  res=goseq(nullP, 'hg38', 'ensGene')
  resUp=res
  resUp$FDR=p.adjust(res$over_represented_pvalue,'fdr')
  resUp=resUp[resUp$FDR<FDR,]
  resUp$Term=sapply(mget(resUp[,1],GOTERM),function(x){x@Term})
  resUp$FoldEnrich=resUp[,4]/resUp[,5]/mean(DE)
  if(!overOnly){
    resDn=res
    resDn$FDR=p.adjust(res$under_represented_pvalue,'fdr')
    resDn=resDn[resDn$FDR<FDR,]
    resDn$Term=sapply(mget(resDn[,1],GOTERM),function(x){x@Term})
    resDn$FoldEnrich=resDn[,4]/resDn[,5]/mean(DE)
    res=rbind(resUp,resDn)
    res$Pvalue=ifelse(res$FoldEnrich>1,res$over_represented_pvalue,res$under_represented_pvalue)
  }else{
    res=resUp
    res$Pvalue=res$over_represented_pvalue
  }
  
  if(addCI){
    if(nrow(res)>0){
      res$nbInGrp=rep(length(geneList),nrow(res))
      res$nbInBckgd=rep(length(background),nrow(res))
      OR=sapply(as.data.frame(t(mapply(function(xy,x,y,tot){odds.ratio(counts2tab(xy,x,y,tot))},res$numDEInCat,res$numInCat,length(geneList),res$nbInBckgd))),as.numeric)
      if(is.matrix(OR)){
        OR=as.data.frame(OR)
      }else{
        nn=names(OR)
        OR=as.data.frame(matrix(OR,1))
        colnames(OR)=nn
      }
      
      OR$CI=paste('[',round(OR[,'LowerCI'],1),'-',round(OR[,'UpperCI'],1),']',sep='')
      res$OR=OR$OR
      res$CI=OR$CI
      res$lowerCI=OR$LowerCI
    }
  }    
  res
}

plotGO=function(resGO,mar=18,nmax=40,colGO=c("#66C2A5","#FC8D62","#8DA0CB"),...){
  names(colGO)=c('BP','MF','CC')
  splitname=function(x,sep=' ',nmax=40){y=strsplit(x,sep)
  y=sapply(y,function(z){
    countchar=cumsum(nchar(z));
    group=countchar%/%nmax;
    paste(by(z,group,paste,collapse=' '),collapse='\n')
  })
  x[nchar(x)>(nmax+10)]=y[nchar(x)>(nmax+10)]
  x
  }														
  params=par()
  par(mar=c(4,mar,1,1))
  resGO=resGO[nrow(resGO):1,]
  barplot(-log10(resGO$FDR), main="", horiz=TRUE, names.arg=splitname(resGO$Term,nmax=nmax),col=ifelse(resGO$ontology=='BP',colGO[1],ifelse(resGO$ontology=='MF',colGO[2],colGO[3])),xlab=expression(-log[10](P[adj])),las=1,...)
  legend("bottomright",fill=colGO[1:3],legend=c('BP','MF','CC'),bty='n')
  par(mar=params$mar)
}