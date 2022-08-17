# NetRankR
This is the GitHub page for the NetRank R package, which provides an implementation for constructing from gene expression data a protein-protein interaction network with StringDB and a weighted gene co-expression network with WGCNA.  


Additionally, this package allows to apply our NetRank algorithm on those networks as described in [this publication](https://www.frontiersin.org/articles/10.3389/fbinf.2022.780229/full).  
Furthermore, additional functionality to analyze clusters within the gene co-expression network, including gene hubs and most common Gene Ontology terms, is provided.  

## Getting started
The R package can simply be installed directly from GitHub by running the following command in R:

`devtools::install_github("Omics-NetRank/NetRankR")`

or by downloading the package as a .zip and using

`devtools::install_local("path/to/package.zip")`

For more details see the provided UserGuide.pdf  

### Citing
If you use the NetRank R package, please cite the following publication:

Al-Fatlawi, A., Afrin, N., Ozen, C., Malekian, N., & Schroeder, M. (2022). NetRank Recovers Known Cancer Hallmark Genes as Universal Biomarker Signature for Cancer Outcome Prediction. Frontiers in Bioinformatics, 2, 780229.

### Contact
  
Ali Al-Fatlawi  
  
Tel: +49 (0)351-463-40063  
Email: Ali.Al-fatlawi@tu-dresden.de  
Webpage: [http://www.biotec.tu-dresden.de/](http://www.biotec.tu-dresden.de/)
