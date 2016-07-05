

DNMFET <- function(clusterslist,runName=""){ 
  if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("parallel needed for this function to work. Please install it.",
    call. = FALSE)
  }
#  library('parallel')

  ### create a matrix for results
  FBET=matrix(nrow=length(clusterslist), ncol=13)
  row.names(FBET)=names(clusterslist)
  colnames(FBET)=c("module size","patho","FET p.value","FET FDR","OR","[95% OR CI] inf","OR [95% OR CI] sup",
                        "module DNMs in patients","module DNMs in controls",
                        "non module DNMs in patients","non module DNMs in controls",
                        "gene names of modules DNMs in patients",
                        "gene names of modules DNMs in controls"
                      )
  
  npDNM=list()
  for(pat in 1:length(PathoGene)){
    cat('\t=====================',names(PathoGene)[pat],'=====================',pat,' of ',length(PathoGene),'\n')
    pathoENS=Patho_ENSgeneID[[pat]]
    DNM_Gene=PathoGene[[pat]]
    ctrlENS=ctr_ENSgeneID[['lgd']]
    ctrl_Gene=ctrGene[['lgd']]
    
    ### function to fill the matrix of results

    #for(i in 1:length(clusterslist)){
    FUNC <- function(i){

      Ms=length(clusterslist[[i]])

      #### FET
      cat('\t\t',names(clusterslist)[i],'\n')
      ## function to calculate the number Mc of DNMs in CTRL involving a gene of the cluster i
      y=lapply(clusterslist[[i]],FUN=function(x) {ctrlENS[which(ctrlENS$ensembl_gene_id==x),'external_gene_name']})
      Mc=sum(sapply(as.matrix(unique(y)), FUN=function(ym) {length(which(ctrl_Gene==ym ))}))
      McID=paste(unlist(y),collapse=", ") 
      
      ## number NMc of remaining DNMs in CTRL involving a gene not in the cluster i
      NMc=length(ctrl_Gene)-Mc
      
      ##function to calculate the number Mee of DNMs in patho involving a gene of the cluster i
      z=lapply(clusterslist[[i]],FUN=function(x) {pathoENS[which(pathoENS$ensembl_gene_id==x),'external_gene_name']})
      Mp=sum(sapply(as.matrix(unique(z)), FUN=function(zm) {length(which(DNM_Gene==zm ))}))
      MpID=paste(unlist(z),collapse=", ") 
      
      ## number NMee of remaining DNMs in EE involving a gene not in the cluster i  
      NMp=length(DNM_Gene)-Mp
      
      # contingency matrice for Fisher Exact Test FET all DNMs and ns DNMs
      matr=matrix(c(Mp,Mc,NMp,NMc), nrow=2)
      
      # FET
      FisherM=fisher.test(matr)
      Fisher.p=FisherM$p.value
      Fisher.or=FisherM$estimate
      Fisher.cinf=FisherM$conf.int[1]
      Fisher.cis=FisherM$conf.int[2]

      FBET[i,]=c(Ms,names(PathoGene[pat]),Fisher.p,NA,Fisher.or,Fisher.cinf,Fisher.cis,Mp,Mc,NMp,NMc,MpID,McID)
    }

    fbet=mclapply(1:length(clusterslist),FUNC,mc.cores=detectCores())
    #The fbet output object of the mclapply function is a list of 44 vectors FBET[i,] in the correct order

    for(i in 1:length(clusterslist)){
      FBET[i,]=fbet[[i]]
    }

    FBET[,"FET FDR"]=p.adjust(FBET[,"FET p.value"],method="fdr")
    npDNM[[pat]]=FBET
  }  
  names(npDNM)=names(PathoGene)

  return(invisible(npDNM))
}


