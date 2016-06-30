



CellTax <- function(clusterslist,outpath=getwd(),runName=""){
  if (!requireNamespace("gplots", quietly = TRUE)) {
    stop("gplots needed for this function to work. Please install it.",
      call. = FALSE)
  }
  #library(gplots)


  pdf(file = paste0(outpath,"/",runName,"_CellTaxonomyHeatmap.pdf"), paper="a4", width = 10,height=15)

  for (i in 1:length(clusterslist)){
    inter=intersect(clusterslist[[i]],rownames(cellexH))
    sel=cellex.an[cellex.an[,"hsapiens_homolog_ensembl_gene"]%in%inter,]
    selH=sel[,c(2:1810)]
    rownames(selH)=sel[,"external_gene_name.y"]
    df=list()
    for(ct in 1:length(ctyp)){
      df[[names(ctyp)[ct]]]=apply(selH[,unlist(ctyp[ct])],1,mean)
    }
    mat=as.matrix(as.data.frame(df))
    mats=t(scale(t(mat)))
    if(anyNA(mats)){
      mats=mats[-which(is.na(mats[,2])==T),]
      }    
    heatmap.2(mats,
              main= names(clusterslist)[i],
              Colv=NULL,
              scale="none",
              dendogram="row",
              col=bluered(200),
              ColSideColors=ctypcol,
              cexRow = 0.15,
              key.title = NA,
              key.xlab = NA,
              key.ylab = NA,
              trace='none',
              #srtCol=0,
              keysize = 1.1,
              density.info='none')
  }
  plot.new()
    legend(x="topright",bty="n",box.lwd = 0.01,box.col = "white",title="GABAergic neurons",
    fill=ctypcol[1:23],
    legend=names(ctyp)[1:23],cex=1)
  plot.new()
    legend(x="topright",bty="n",box.lwd = 0.01,box.col = "white",title="Glutamatergic neurons",
    fill=ctypcol[24:42],
    legend=names(ctyp)[24:42],cex=1)
  plot.new()
    legend(x="topright",bty="n",box.lwd = 0.01,box.col = "white",title="Non neuronal cells",
    fill=ctypcol[43:49],
    legend=names(ctyp)[43:49],cex=1)
  
  dev.off()

}
