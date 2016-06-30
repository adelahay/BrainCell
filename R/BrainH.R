BrainH <- function(clusterslist,outpath=getwd(),runName=""){  
  
  NfuGSE25219=rbind(NfuGSE25219_1,NfuGSE25219_2)
  
  if (!requireNamespace("gplots", quietly = TRUE)) {
  stop("gplots needed for this function to work. Please install it.",
    call. = FALSE)
  }
  #library(gplots)
  #load(paste0(inpath,'/NfuGSE25219ensg.Rdata')) #
  #load(paste0(inpath,'/GSE25219forHeatmap_GSM.Rdata')) #GSM


  pdf(file = paste(outpath,"/",runName,"_BainHeatmap.pdf",sep=""), paper="a4r", width = 10)
  for (i in 1:length(clusterslist)){
    NAmat <- matrix(data=as.numeric(NA),nrow=nrow(GSM),ncol=ncol(GSM),dimnames= dimnames(GSM))
      ME.av = function(mat,data){
        for (r in 1:nrow(mat)){
          for(c in 1:ncol(mat)){
            mat[r,c] <- mean(data[unlist(GSM[r,c])])
          }
        }
        return(mat)
      }
    
    Mgenes= intersect(clusterslist[[i]],rownames(NfuGSE25219))
    
    if(length(Mgenes)>=1){
      if(length(Mgenes)==1){
        Mdata <- NfuGSE25219[Mgenes,]
      }else{
        Mdata <- apply(NfuGSE25219[Mgenes,],2,mean) 
      }
      Mavmat <- ME.av(NAmat,Mdata) 
      heatmap.2(Mavmat,
                main= names(clusterslist)[i],
                Colv=NULL,
                col=bluered(200),
                key.title = NA,
                key.xlab = NA,
                key.ylab = NA,
                trace='none',
                srtCol=0,
                keysize = 1.1,
                density.info='none')
    }
  }
  dev.off()
  

}




