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

