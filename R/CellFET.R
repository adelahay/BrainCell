CellFET <- function(clusterslist,runName=""){	#inpath="~/Dropbox/SHARED/tools/Data_to_load_CellFET",

#cat("Cell background is the genes with one2one human orthologous of mice genes used to build the list of cell class enriched genes by Zeisel et al 2015(Science)")

if (!requireNamespace("parallel", quietly = TRUE)) {
stop("parallel needed for this function to work. Please install it.",
  call. = FALSE)
}
#  library('parallel')
  
  ### create a matrix for results
  cFET=matrix(nrow=length(clusterslist), ncol=14)
  row.names(cFET)=names(clusterslist)
  colnames(cFET)=c("cell class","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                      "module cell enriched genes","module out cell enriched genes",
                      "non module cell enriched genes","non module out cell enriched genes",
                      "gene names of in modules cell enriched genes","module size","cell enriched genes size",
                      "% of genes in module in cell background"
  )

  
  resMsc=list()
  for(ccl in 1:length(hmscDF)){ # ccl: cell class
    cat('\t=====================',names(hmscDF)[ccl],'=====================',ccl,' of ',length(hmscDF),'\n')
    cclENS=hmscDF[[ccl]]
    ### function to fill the matrix of results
    #for(i in 1:length(clusterslist)){
    FUNC=function(i){
      Ms=length(clusterslist[[i]]) #Ms: module size
      CB=HUMQb[,'hsapiens_homolog_ensembl_gene'] #CB: cell background
      Cs=length(cclENS) #Cs: cell enriched genes size
      MCBp=length(intersect(CB,clusterslist[[i]]))/Ms #MCBp: % of genes in module in cell background

      #cFET
      cat('\t\t',names(clusterslist)[i],'\n')
      #calculate the number Mc of module i cell enriched genes(Mc: in module AND in cell class)
      Mc=length(intersect(cclENS,clusterslist[[i]]))
      McID=paste(unlist(HUMQb[which(CB %in% intersect(cclENS,clusterslist[[i]])),'external_gene_name']),collapse=", ")
      #calculate the number NMc of remaining genes not in module but in cell class
      NMc=length(cclENS)-Mc
      #calculate the number Mnc of genes in module but not in cell class
      Mnc=length(intersect(CB,clusterslist[[i]]))-Mc
      #calculate the number NMnc of genes out of module AND not in cell class
      NMnc=length(CB)-(Mc+NMc+Mnc)
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr=matrix(c(Mc,NMc,Mnc,NMnc), nrow=2)
      #FET
      #FisherM=fisher.test(matr,alternative="greater")
      FisherM=fisher.test(matr)
      Fisher.p=FisherM$p.value
      Fisher.or=FisherM$estimate
      Fisher.cinf=FisherM$conf.int[1]
      Fisher.cis=FisherM$conf.int[2]
      cFET[i,]=c(names(hmscDF)[ccl],Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mc,Mnc,NMc,NMnc,McID,Ms,Cs,MCBp)
    }
#    cfet=mclapply(1:length(clusterslist),FUNC,mc.cores=detectCores())
    cfet=lapply(1:length(clusterslist),FUNC)
    #The cfet output object of the mclapply function is a list of n vectors cFET[i,] in the correct order
    for(i in 1:length(clusterslist)){
      cFET[i,]=cfet[[i]]
    }
    cFET[,"FET FDR"]=p.adjust(cFET[,"FET p.value"],method="fdr")
    resMsc[[ccl]]=cFET
  }  
  names(resMsc)=names(hmscDF)

  return(resMsc)      
}

