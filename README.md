# NetRankR
This is the GitHub page for the NetRank R package, which provides an implementation for constructing from gene expression data a protein-protein interaction network with StringDB and a weighted gene co-expression network with WGCNA.  


Additionally, this package allows to apply our NetRank algorithm on those networks as described in [this publication](https://www.frontiersin.org/articles/10.3389/fbinf.2022.780229/full).  
Furthermore, additional functionality to analyze clusters within the gene co-expression network, including gene hubs and most common Gene Ontology terms, is provided.  

## Getting started
The R package can simply be installed directly from GitHub by running the following command in R:

`devtools::install_github("Omics-NetRank/NetRankR")`

or by downloading the package as a .zip and using

`devtools::install_local("path/to/package.zip")`

For more details see the provided [UserGuide.pdf](https://github.com/Alfatlawi/Omics-NetRank/blob/master/UserGuide.pdf).


### Test data
We evaluate our algorithm and implementation on data for 3,388 patients for 19 cancer types gathered from the Cancer Genome Atlas (TCGA) (https://portal.gdc.cancer.gov/). 
We obtained gene expression data from the Cancer Genome Atlas (TCGA), which initially consisted of 20,531 genes and 11,069 samples. We kept only 8,603 samples after removing samples with missing values and duplicates. Among these, we used only 3,388 samples covering 19 cancer types, which were manually reviewed and approved in their clinical follow-up (see data distribution plot in the images folder ).  We normalized the expression data by using the MinMaxScaler function of the scikit-learn package version 1.0.2 (Pedregosa, F., Varoquaux, G., Gramfor...). For each cancer type, we split the data into 70/30\% sets, using the 70\% for network construction and gene ranking generation and the remaining 30\% for PCA and SVM testing and evaluation.

For the phenotype and the normalized expression data, please refer to this [link]([https://cloudstore.zih.tu-dresden.de/index.php/s/77xCWyqnpStFPjL](https://sharing.biotec.tu-dresden.de/index.php/s/Fw7sjbbt9jfZWGx).

...

### Contact
  
Ali Al-Fatlawi  
  
Tel: +49 (0)351-463-40063  
Email: Ali.Al-fatlawi@tu-dresden.de  
Webpage: [http://www.biotec.tu-dresden.de/](http://www.biotec.tu-dresden.de/)
