require(ggplot2)
require(NMF)
require(RColorBrewer)

#Create scale for heatmap
hmpscaletest = function(hmpcol,voom,checkbox){
  min=min(voom)
  max=max(voom)
  val=sort(c(seq(min,max, length.out = 5)),decreasing=TRUE)
  df <- data.frame(x = rep(1, 5),y = val,z = factor(1:5))
  if(checkbox==FALSE){
    g1=ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
  }
  else{
    g1=ggplot(df, aes(x, y)) +geom_tile(aes(fill = z))+scale_fill_brewer( type = "div" , palette = hmpcol,direction = -1)+guides(fill=FALSE)+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+coord_flip()
  }
  return(g1)
}
###########################################################
#Function to generate limma heatmaps's
createheatmap <- function(results,expr,pval,hmpsamp,contrast){
  if(missing(pval)){
    top_expr=expr
  }else{
top_expr=expr[match(rownames(pval),rownames(expr)),]}
if(hmpsamp==F){
  contr=strsplit(contrast,"_vs_")
  ct1=sapply(contr,"[",1)
  ct2=sapply(contr,"[",2)
  pd=pData(results$eset)
  sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
  sample=as.character(sample)
  top_expr=top_expr[,eval(sample)]}
return(top_expr)
}

###########################################################

#Function to generate camera, spia and go heatmap
heatmapfun <- function(results,expr,pval,file,prj,hmplim,hmpsamp,contrast){
  dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}
  pd=pData(results$eset)
  old=file$old[file$projects %in% prj]
  old=as.character(old)
  seq=file$seq[file$projects %in% prj]
  seq=as.character(seq)
  if(seq=="R"){
    expr$ENSEMBL=rownames(expr)
    expr=inner_join(expr,pval,by=c('ENSEMBL'='ENSEMBL'))
    expr=expr[order(expr$adj.P.Val),]
    rownames(expr)=make.names(expr$SYMBOL,unique=T)
    if(old=="N"){
      expr=expr %>% dplyr::select(-ENSEMBL:-t)}
    else if(old=="Y"){
      expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
    }
  }
  else if(seq=="M"){
    expr$id=rownames(expr)
    pval$id=rownames(pval)
    expr=inner_join(expr,pval,by=c('id'='id'))
    expr=expr[order(expr$adj.P.Val),]
    rownames(expr)=make.names(expr$SYMBOL,unique=T)
    if(old=="N"){
      expr=expr %>% dplyr::select(-ENSEMBL:-t)
    }
    else if(old=="Y"){
      expr=expr %>% dplyr::select(-ENSEMBL:-adj.P.Val)
    }
  }
  
  top_expr=as.data.frame(expr)
  top_expr=top_expr[1:hmplim,]
  samples=as.character(rownames(pd))
  top_expr=top_expr %>% dplyr::select(samples)
  if(hmpsamp==F){
    contr=strsplit(contrast,"_vs_")
    ct1=sapply(contr,"[",1)
    ct2=sapply(contr,"[",2)
    sample=pd$sample_name[pd$maineffect %in% c(ct1,ct2)]
    sample=as.character(sample)
    top_expr=top_expr[,eval(sample)]}
  return(top_expr)
}

###########################################################
#functions to create heatmaps and other helper functions
library(DOSE)
library(stats)
library(qvalue)
library(fgsea)
library(enrichplot)
creategseaobj = function(geneList, geneSets, exponent=1,
                         nPerm=1000,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff=1,
                         pAdjustMethod="BH"){
  
  tmp_res <- fgsea::fgsea(pathways=geneSets,
                          stats=geneList,
                          nperm=nPerm,
                          minSize=minGSSize,
                          maxSize=maxGSSize,
                          gseaParam=exponent,
                          nproc = 0)
  
  p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
  qvalues <- calculate_qvalue(tmp_res$pval)
  
  params <- list(pvalueCutoff = pvalueCutoff,
                 nPerm = nPerm,
                 pAdjustMethod = pAdjustMethod,
                 exponent = exponent,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize
  )
  res <- data.frame(
    ID = as.character(tmp_res$pathway),
    Description = tmp_res$pathway,
    setSize = tmp_res$size,
    enrichmentScore = tmp_res$ES,
    NES = tmp_res$NES,
    pvalue = tmp_res$pval,
    p.adjust = p.adj,
    qvalues = qvalues,
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$pvalue),]
  res <- res[ res$pvalue <= pvalueCutoff, ]
  res <- res[ res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  
  row.names(res) <- res$ID
  observed_info <- lapply(geneSets[res$ID], function(gs)
    gseaScores(geneSet=gs,
               geneList=geneList,
               exponent=exponent)
  )
  
  #   if (verbose)
  #     message("leading edge analysis...")
  
  ledge <- leading_edge(observed_info)
  
  res$rank <- ledge$rank
  res$leading_edge <- ledge$leading_edge
  res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')
  
  #   if (verbose)
  #     message("done...")
  
  new_res=new("gseaResult",
              result     = res,
              geneSets   = geneSets,
              geneList   = geneList,
              params     = params,
              readable   = FALSE
  )
  
  return(new_res)
}

###########################################################
calculate_qvalue <- function(pvals) {
  if (length(pvals) == 0)
    return(numeric(0))
  
  qobj <- tryCatch(qvalue(pvals, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)
  
  if (class(qobj) == "qvalue") {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }
  return(qvalues)
}

###########################################################
gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}
###########################################################

leading_edge <- function(observed_info) {
  core_enrichment <- lapply(observed_info, function(x) {
    runningES <- x$runningES
    runningES <- runningES[runningES$position == 1,]
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      leading_gene <- runningES$gene[1:i]
    } else {
      i <- which.min(runningES$runningScore)
      leading_gene <- runningES$gene[-c(1:i)]
    }
    return(leading_gene)
  })
  rank <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    ES <- x$ES
    if (ES >= 0) {
      rr <- which.max(runningES$runningScore)
    } else {
      i <- which.min(runningES$runningScore)
      rr <- nrow(runningES) - i + 1
    }
    return(rr)
  })
  
  tags <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    runningES <- runningES[runningES$position == 1,]
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      res <- i/nrow(runningES)
    } else {
      i <- which.min(runningES$runningScore)
      res <- (nrow(runningES) - i + 1)/nrow(runningES)
    }
    return(res)
  })
  
  ll <- sapply(observed_info, function(x) {
    runningES <- x$runningES
    ES <- x$ES
    if (ES >= 0) {
      i <- which.max(runningES$runningScore)
      res <- i/nrow(runningES)
    } else {
      i <- which.min(runningES$runningScore)
      res <- (nrow(runningES) - i + 1)/nrow(runningES)
    }
    return(res)
  })
  
  N <- nrow(observed_info[[1]]$runningES)
  setSize <- sapply(observed_info, function(x) sum(x$runningES$position))
  signal <- tags * (1-ll) * (N / (N - setSize))
  
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <- paste0('tags=', tags, ", list=", ll, ", signal=", signal)
  
  res <- list(rank = rank,
              tags = tags,
              list = ll,
              signal = signal,
              leading_edge = leading_edge,
              core_enrichment = core_enrichment)
  return(res)
}

###########################################################
#findinf right geneset for the data
findgeneset <- function(organism,geneset){
  hum=c("human","Human","Hs","Homo sapiens","Homo_sapiens")
  mouse=c("mouse","Mouse","Mm","Mus musculus","Mus_musculus")
  rat=c("rat","Rat","Rn","Rattus norvegicus","Rattus_norvegicus")
  if(organism %in% hum){
    if(geneset =="Hallmark"){
    data('human_H_v5',package="ExpressExtras")
    geneset = Hs.H
    }else if(geneset =="Curated"){
    data('human_c2_v5',package="ExpressExtras")
    geneset = Hs.c2
    }
  }else if(organism %in% mouse){
    if(geneset =="Hallmark"){
    data('mouse_H_v5',package="ExpressExtras")
      geneset = Mm.H
    }else if(geneset =="Curated"){
    data('mouse_c2_v5',package="ExpressExtras")
      geneset = Mm.c2
    }else if(geneset =="GO"){
    data('mouse_GO',package="ExpressExtras")
      geneset = Mm.GO
    }
  }else if(organism %in% rat){
    if(geneset =="Hallmark"){
    data('Rat_Hallmark',package="ExpressExtras")
      geneset = Rn.H
    }else if(geneset =="Curated"){
    data('Rat_C2_v4p2',package="ExpressExtras")
      geneset = Rn.c2
    }else if(geneset =="GO"){
    data('Rat_GO',package="ExpressExtras")
      geneset = Rn.GO
  }
    }else {
    print("Incorrect organism name. ")}
  return(geneset)
}