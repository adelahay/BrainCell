#### BrainCell

The BrainCell R package provides functions to test query gene sets in the context of brain and neurodevelopment
     
  <i>BrainH</i> represents gene sets expression pattern across brain regions during development and post-natal life using the [GSE25219 dataset](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25219). The acronyms for brain regions and the numeration of lifetime periods were previously designated by [Kang et al.](http://www.ncbi.nlm.nih.gov/pubmed/22031440) and are used to label the heatmaps.  
     
  <i>CellTax</i> represents gene sets expression pattern across cortical cell types using the [GSE71585 dataset](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585). The cell type nomenclature were designated by [Tasic et al.](http://www.ncbi.nlm.nih.gov/pubmed/26727548) and are used to label the heatmaps.   
     
  <i>CellFET</i> tests gene sets for enrichment in brain cell type marker genes for Ependymal cells, Oligodendrocytes, Microglial cells, Pyramidal neurons from CA1 area, Interneuron, Endothelial cells, Pyramidal neurons from S1 area, Astrocytes and Mural cells (as designated by [Zeisel et al.](http://www.ncbi.nlm.nih.gov/pubmed/25700174))  
     
  <i>DNMFET</i> tests gene sets for enrichment in deleterious <i>de novo</i> mutations ascertained from patients with neurodevelopmental disorder (EE= Epileptic Encephalopathy, ASD= Autism Spectrum Disorder, ID=Intellectual Disability, SCZ= Schizophrenia, combined= EE+ASD+ID+SCZ, DDD= Developmental disorder from the DDD study) as previously described by [Johnson et al.](http://www.ncbi.nlm.nih.gov/pubmed/26691832)  


#### Installation

You can install the BrainCell package from github.

    install.packages('devtools')
    library(devtools)

    install_github("adelahay/BrainCell")
    library(BrainCell)
    
In addition to this, you would need to install all package dependencies. They are downloadable from CRAN (http://cran.r-project.org)

    install.packages("gplots")  
    library('gplots')
    install.packages("parallel")
    library('parallel')
    


#### Examples/Tutorial

The query gene sets need to be a list of character vectors of human ENSEMBL gene ID. In the example we run the function for <i>CLtest</i>, a list of two gene sets: the first gene set "M30" contains 320 genes and the second "M18" contains 149 genes. 

     is(CLtest)
     #[1] "list"   "vector"
     names(CLtest)
     #[1] "M30" "M18"
     length(CLtest[[1]])
     #[1] 320
     head(CLtest[[1]])
     #"ENSG00000183527" "ENSG00000145526" "ENSG00000119946" "ENSG00000151849" "ENSG00000184524" "ENSG00000081148"
     length(CLtest[[2]])
     #[1] 149
     
If we create a new folder "/BrainCellH_test_outputs" for outputs and save its path as "testpath"

     system(paste0("mkdir ",getwd(),"/BrainCellH_test_outputs"))
     testpath=paste0(getwd(),"/BrainCellH_test_outputs")

We run the BrainCell functions

     BrainH(clusterslist = CLtest, outpath = testpath, runName = "test")

<i>BrainH</i> function produces a pdf [test_BrainHeatmap.pdf](https://github.com/adelahay/BrainCell/blob/master/BrainCellH_test_outputs/test_BainHeatmap.pdf)


     CellTax(clusterslist = CLtest, outpath = testpath, runName = "test")

<i>CellTax</i> function produces a pdf [test_CellTaxonomyHeatmap.pdf](https://github.com/adelahay/BrainCell/blob/master/BrainCellH_test_outputs/test_CellTaxonomyHeatmap.pdf) 
     

     resC=CellFET(clusterslist = CLtest, runName = "test")
     resD=DNMFET(clusterslist = CLtest, runName = "test")

<i>CellFET</i> and <i>DNMFET</i> functions need to be adressed to a results R object (here, resC and resD). These results R objects are a list of matrices with one matrice for each cell type tested in <i>CellFET</i> function, and for each phenotype tested in <i>DNMFET</i> function.  
   

The following code will produce .txt files to have a overview of the results: [test_cellFET.txt](https://github.com/adelahay/BrainCell/blob/master/BrainCellH_test_outputs/test_cellFET.txt) and [test_DNMFET.txt](https://github.com/adelahay/BrainCell/blob/master/BrainCellH_test_outputs/test_DNMFET.txt)  

     # to save the CellFET outputs in .txt file
     allT<- resC[[1]]
          for (ccl in 2:length(resC)){
               allT<- rbind(allT,resC[[ccl]])
          }
     write.table(allT, sep = '\t', file =paste0(testpath,'/test_cellFET.txt'), row.names = TRUE, quote = FALSE, col.names = NA)
     

     # to save the DNMFET outputs in .txt file
     allT<- resD[[1]]
          for (ccl in 2:length(resD)){
               allT<- rbind(allT,resD[[ccl]])
          }
     write.table(allT, sep = '\t', file =paste0(testpath,'/test_DNMFET.txt'), row.names = TRUE, quote = FALSE, col.names = NA)


#### Citation

Please cite [Delahaye et al. 2016 Genome Biology (link to complete)](http://dx.doi.org/)

[![DOI](https://zenodo.org/badge/62316752.svg)](https://zenodo.org/badge/latestdoi/62316752)


