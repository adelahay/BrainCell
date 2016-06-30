#### BrainCell

The BrainCell R package provides functions to test query gene sets in the context of brain and neurodevelopment
     
  <i>BrainH</i>  
  represents gene sets expression pattern across brain regions during development and post-natal life  
     
  <i>CellTax</i>  
  represents gene sets expression pattern across cortical cell types  
     
  <i>CellFET</i>  
  tests gene sets for enrichment in brain cell type marker genes  
     
  <i>DNMFET</i>  
  tests gene sets for enrichment in deleterious <i>de novo</i> mutations ascertained from patients with neurodevelopmental disorders  
     
  <i>forest.plot</i>  
  plots the results of <i>CellFET</i> and <i>DNMFET</i> functions


#### Installation

You can install the BrainCell package from github.

    install.packages('devtools')
    library(devtools)

    install_github("adelahay/BrainCell")
    library(BrainCell)
    
In addition to this, you would need to install all package dependencies. They are downloadable from CRAN (http://cran.r-project.org).

    #Installing from CRAN
    install.packages("gplots")  
    library('gplots')
    install.packages("ggplot2")
    library('ggplot2')
    install.packages("parallel")
    library('parallel')
    


#### Examples/Tutorial

The query gene sets need to be a list of character vectors of human ENSEMBL gene ID. In the example we run the function for a list of two gene sets: the first gene set "M30" contains 320 genes and the second "M18" contains 149 genes. 

#### Citation

Please cite [link to complete](http://dx.doi.org/)


