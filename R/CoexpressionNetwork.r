##
#' Outlier Detection for WGCNA
#' @description Creates files to manually detect outliers for creating a the coexpression network. Plots are saved to workingDir/datasetDir/Plots.
#' @param networkType Type of coexpression network to be created. Default = "unsinged".
#' @param datasetDir Name of the directory containing the Dataset. The Directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param pathToExprData Path to the csv file which contains the expression data.
#' @param pathToTraitData Path to the csv file which contains the phenotype data.
#' @param sampleIdCol Column name of sample IDs. Default = "samplename".
#' @param phenoCol Column name of phenotypes. Default = "Pheno".
#' @return Creates and saves pdf files in workingDir/dataset/networkType/Plots/
#' @examples
#' check README for an example
#' @export
CoexpressionNetwork_Preprocessing <- function(networkType = "unsigned", datasetDir, workingDir, pathToExprData, pathToTraitData, sampleIDCol = 'sampleName', phenoCol = 'pheno') {


  setwd(workingDir);
  doParallel::registerDoParallel(cores=30)
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);

  dir.create(file.path(workingDir, datasetDir, networkType), showWarnings = FALSE)
  dir.create(file.path(workingDir, datasetDir, networkType, "Plots"), showWarnings = FALSE)

  # read data
  datExpr0 <- read.csv(pathToExprData, sep=",");
  datTraits <- read.csv(pathToTraitData, sep=",");
  rownames(datExpr0) <- datExpr0[,1]
  datExpr0 = as.data.frame((datExpr0[,-c(1)]),check.names=FALSE); #


  #remove duplicates
  datExpr0 = as.data.frame(t(datExpr0))
  n_occur <- data.frame(table(rownames(datExpr0)))

  duplicatedGenes = n_occur[n_occur$Freq > 1,]
  genes = duplicatedGenes[,1]
  freq = duplicatedGenes[,1]
  index = rownames(duplicatedGenes)
  datExpr0 = subset(datExpr0, !(rownames(datExpr0) %in% index))

  # remove columns that hold information we do not need.
  allTraits = datTraits[, c(sampleIDCol, phenoCol)]
  names(allTraits) = c("sampleName", "pheno")
  rownames(allTraits) = allTraits$sampleName
  allTraits = allTraits[allTraits$pheno!=-9,]

  ## keep only the samples that have data and labels:
  intersectedSamples = intersect(colnames(datExpr0), rownames(allTraits))
  datTraits0 = allTraits[intersectedSamples,]
  datExpr0 = datExpr0[,intersectedSamples]

  # Form a data frame analogous to expression data that will hold the clinical traits. Match the order of samples and phenotype
  Samples = names(datExpr0);
  traitRows = match(Samples, datTraits0$sampleName);
  traitRows = traitRows[!is.na(traitRows)]
  allTraits = allTraits[traitRows,]

  # convert phenotype into numerical:
  datTraits0$pheno = factor(datTraits0$pheno)
  # Iterative Filtering Of Samples And Genes With Too Many Missing Entries
  gsg = WGCNA::goodSamplesGenes(t(datExpr0), verbose = 3);
  gsg$allOK

  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      print(paste("Removing genes:", paste(length(colnames(t(datExpr0))[!gsg$goodGenes]), collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      print(paste("Removing samples:", paste(length(rownames(t(datExpr0))[!gsg$goodSamples]), collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodGenes, gsg$goodSamples]

    # Maybe some samples were discarded in the previous step
    Samples = colnames(datExpr0);
    traitRows = match(Samples, datTraits0$sampleName);
    traitRows= traitRows[!is.na(traitRows)]
    s0=datTraits0[traitRows,]
    dim(datTraits0)
    dim(datExpr0)
    print('dim after remove bad genes or samples')
    print(dim(datExpr0))
  }



  if (all(colnames(datExpr0) == datTraits0[, 1])){print('Message: the traits and expression data have been aligned correctly.')
  }else {print("Error Message: the traits and expression data are not aligned correctly. This should not happen, please check your data.")
  }


  # transpose and remove annotation:
  #We now remove the auxiliary data and transpose the expression data for further analysis.
  datExpr0=as.data.frame(t(datExpr0),check.names=FALSE);

  #Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
  print("Average clustering in progress")
  sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  # The plots are saved within workingDir/datasetDir/networkType/Plots
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_sampleClustering.pdf",sep = ""), width = 12, height = 9);par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_sampleClustering_pheno.pdf",sep = ""), width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", labels=datTraits0$pheno, sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()


  #### PCA
  # Compute PCA
  print("Computing PCA")
  res.pca <- prcomp(datExpr0, scale = F)

  #Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_PCA_eigenvalues.pdf",sep = ""), width = 12, height = 9);
  print(factoextra::fviz_eig(res.pca))
  dev.off()

  #Graph of individuals. Individuals with a similar profile are grouped together.
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_PCA_individuals.pdf",sep = ""), width = 12, height = 9);

  print(factoextra::fviz_pca_ind(res.pca,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
  ))
  dev.off()
  ## Plots to help choose soft thresholding power
   #Figure 1: Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit
  #index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity
  #(degree, y-axis) as a function of the soft-thresholding power (x-axis).

  print("Saving now!")

  save(datExpr0, datTraits0, sampleTree, file= paste(datasetDir, "/CoexpressionNetwork_Preprocessing.RData",sep = ""))
  print(paste('First part is done now! Check the plots found in ',workingDir, datasetDir, '/', networkType, "/Plots to decide cut height and power threshold.", sep="" ))
  print('More information regarding cut height and power thresholding can be found in the WGCNA documentation.')

}


##
#' Building the coexpression network
#' @description Creates and stores a coexpression network from gene expression data based on WGCNA blockwiseModules.
#' @param networkType Type of Network to be created. Default = "unsinged".
#' @param power_thre Power Threshold for soft thresholding during network construction. The current version of NetRank does not include functionality to help the user choose the power threshold. For more information see the WGCNA publication.
#' @param maxBlockSize Maximum block size. Higher values demand more RAM.
#' @param cutHeight Cut height for culling outliers. Choose value by checking plots created by CoexpressionNetwork_Preprocessing. For further information see the WGCNA publication
#' @param datasetDir Name of the directory containing the Dataset. The Directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param cores Number of cores to be used for parallel computation
#' @return Creates a coexpression network and saves it tn the "datasetDir" directory.
#' @examples
#' check README for an example
#' @export
CoexpressionNetwork_BuildNetwork <- function(networkType = "unsigned", power_thre, maxBlockSize, cutHeight, workingDir, datasetDir, cores){

  setwd(workingDir)
  doParallel::registerDoParallel(cores=cores)
  if (cores > 1) {WGCNA::enableWGCNAThreads(nThreads = cores)}

  options(stringsAsFactors = FALSE);

  loadData <- load(file = paste(workingDir, datasetDir,'/CoexpressionNetwork_Preprocessing.RData',sep = ""))

  # Plot a line to show the cut.
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_sampleClustering_lined.pdf",sep = ""), width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  abline(h = cutHeight, col = "red");
  dev.off()

  # Determine cluster under the line
  clust = WGCNA::cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
  table(clust)
  # clust 1 contains the samples we want to keep.
  keepSamples = (clust==1)
  datExpr = datExpr0[keepSamples, ]
  datTraits =datTraits0[keepSamples,]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  #The variable datExpr now contains the expression data ready for network analysis.

  #Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
  sampleTree = hclust(dist(datExpr), method = "average");

  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.

  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_sampleClustering_afterRemove.pdf",sep = ""), width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_sampleClustering_pheno_afterRemove.pdf",sep = ""), width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", labels=datTraits$pheno, sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()

  WGCNA::collectGarbage();

  #We now have the expression data in the variable datExpr, and the corresponding clinical traits in the variable datTraits.
  #Before we continue with network construction and module detection,
  #we visualize how the clinical traits relate to the sample dendrogram.
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = WGCNA::numbers2colors((as.numeric(datTraits$pheno)), signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/1_Sample_dendrogram_trait_heatmap.pdf",sep = ""), width = 12, height = 9);
  WGCNA::plotDendroAndColors(sampleTree2, traitColors,  groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
  dev.off()
  #In the plot, shown in Fig. 2, white means a low value, red a high value, and grey a missing entry.

  #The last step is to save the relevant expression and trait data for use in the next steps

  #If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
  # Iterative Filtering Of Samples And Genes With Too Many Missing Entries
  gsg = WGCNA::goodSamplesGenes(datExpr, verbose = 3);
  gsg$allOK

  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      print(paste("Removing genes:", paste(length(names(datExpr)[!gsg$goodGenes]), collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      print(paste("Removing samples:", paste(length(rownames(datExpr)[!gsg$goodSamples]), collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }



  # Maybe some samples were discarded in the previous step
  Samples = rownames(datExpr);
  traitRows = match(Samples, rownames(datTraits));
  traitRows = traitRows[!is.na(traitRows)]
  datTraits = datTraits[traitRows,]
  dim(datTraits)
  dim(datExpr)

  if (all(rownames(datExpr) == datTraits[, 1])){print('Message: the traits and expression data have been aligned correctly.')
  }else {print("Error Message: the traits and expression data are not aligned correctly. This should not happen, please check your data.")}


  cor <- WGCNA::cor
  TOMType = networkType
  setwd(paste(workingDir, datasetDir, sep='' ))
  if (networkType=='signed hybrid') {TOMType='signed'}

  net = WGCNA::blockwiseModules(datExpr, power = power_thre, networkType=networkType,
                         TOMType = TOMType, minModuleSize = min(20, ncol(datExpr)/2),
                         reassignThreshold = 1e-6, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "CoexpressionNetwork_TOM",maxBlockSize = maxBlockSize,
                         verbose = 3)

  save(net, datExpr, datTraits,  file= paste(networkType, "/CoexpressionNetwork.RData",sep = ""))

  print(paste("Step finished. Results have been saved to ", workingDir, "/", datasetDir, "/CoexpressionNetwork_TOM-block.1.Rdata and ", workingDir, datasetDir,'/', networkType, "/CoexpressionNetwork.RData", sep=""))
  print("If more than one Tom-block file has been created, increase maxBlockSize and rerun this step.")

}

#' Preparing Coexpression Network for NetRank.
#' @description Creates files for NetRank from the coexpression network
#' @param datasetDir Name of the directory containing the dataset. The directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param networkType Type of the previously computed coexpression network, default = "unsigned".
#' @param standardCorr_column Name of the column containing the correlation information in the correlation file. Default = 'PearsonCorrelation'.
#' @param geneName_col Name of the column containing the gene names in the correlation file. Default = 'GeneName'
#' @param dataFolder Name of the directory containing the correlation file. The NetRank results will be saved into this directory. Default = 'standardGeneScreeningResults'.
#' @param NetWorkCutThreshold In the coexpression network edges with a coexpression value below the threshold are discarded. This is necessary due to the size of the network, without the cut it is too large for further processing.
#' @param corrFilename Name of the file containing the gene correlation information.
#' @return Saves an .RData file in the datasetDir directory for NetRank application.
#' @examples
#' check the README for an example
#' @export
CoexpressionNetwork_Postprocessing <- function(datasetDir, workingDir, networkType = 'unsigned', standardCorr_column = 'PearsonCorrelation', geneName_col = "GeneName", dataFolder = 'standardGeneScreeningResults', NetWorkCutThreshold, corrFilename){


  options(stringsAsFactors = FALSE);

  setwd(paste(workingDir, datasetDir, sep='' ))
  print("Loading network information")
  loadWGCNA <- load(file = paste(workingDir,'/', datasetDir, '/CoexpressionNetwork_TOM-block.1.RData', sep = "")) # load WGCNA Network info
  #load data from previous steps
  loaddata <- load(file = paste(workingDir, datasetDir, '/',networkType,"/CoexpressionNetwork.RData",sep = ""))
  # Read the standard_corr from your correlation analysis.
  standard_corr = read.csv(paste(workingDir, datasetDir, "/", dataFolder, '/',corrFilename, sep = ""), sep=",");
  standard_corr <- standard_corr[complete.cases(standard_corr[, standardCorr_column]),] # remove rows with missing or Na value in standardCorr_column
  standard_corr[standardCorr_column] = abs(standard_corr[standardCorr_column]) # Sort by correlation
  standard_corr = standard_corr[c(geneName_col, standardCorr_column)] ## cull nonessential columns
  colnames(standard_corr) = c("GeneName", standardCorr_column)

  #convert network information from the coexpression network construction into data.frame
  adjmatrix = as.matrix(TOM)
  rownames(adjmatrix) = c(colnames(datExpr))
  colnames(adjmatrix) = c(colnames(datExpr))
  tmp <- t(combn(colnames(adjmatrix), 2)) # only upper triangle
  adjmatrix <- data.frame(tmp, corr=adjmatrix[tmp]) # convert matrix to data.frame
  rm(tmp)

  colnames(adjmatrix)=c('from','to', 'corr')

  # reformat the data to fit the NetRank algorithm.
  print("preparing Data for Netrank, for large datasets this will take some time.")
  get_interactions2 = adjmatrix
  rm(adjmatrix)
  rm(datExpr)
  rm(datTraits)
  rm(TOM)
  rm(net)
  rm(loadWGCNA)
  rm(loaddata)
  
  get_interactions2 <- subset(get_interactions2, get_interactions2[ ,3] > NetWorkCutThreshold) # remove weakly correlated edges to reduce size
  get_interactions2 = subset(get_interactions2, is.element(get_interactions2$from, standard_corr[[geneName_col]])) # remove genes without correlation value
  get_interactions2 = subset(get_interactions2, is.element(get_interactions2$to, standard_corr[[geneName_col]]))
  get_interactions2 = get_interactions2[!duplicated(cbind(get_interactions2$from,get_interactions2$to)),  ] # remove duplicated
  get_interactions2 = get_interactions2[!duplicated(cbind(get_interactions2$to,get_interactions2$from)),  ] # remove duplicated
  write.table(get_interactions2, file = paste(workingDir, "/", datasetDir, "/getinteractions.csv",sep = ""), row.names = F)

  ## Take the uniques ids and their counts
  #Ids = unique(get_interactions2$to)
  Ids = unique(c(cbind(get_interactions2$from, get_interactions2$to))) ## get the id's of nodes that connect to the node i.
  N = length(Ids)

  print("Aggregating connectivity scores")
  ### sum of connectivity (outputs) of each genes.
  outs <- data.frame(matrix(ncol = 0, nrow = N))
  outs$from = Ids
  gc(verbose = F)
  outsCount = aggregate(corr ~ from, get_interactions2, sum) # sum of connectivity
  outsCount2 = aggregate(corr ~ to, get_interactions2, sum) # sum of connectivity
  outsCount= merge(outsCount, outsCount2,by.x = 'from', by.y = 'to' , all=T)
  outsCount[is.na(outsCount)] <- 0
  outsCount$sum=outsCount$corr.x +outsCount$corr.y
  outsCount=outsCount[,c('from','sum')]
  outs = merge(outs, outsCount, all.x=T)


  colnames(outs) = c(geneName_col,'outs_conn')
  ## Merge the sum as a column
  standard_corr = merge(standard_corr, outs, by = geneName_col)

  rm(outs)
  rm(outsCount)

  save(standard_corr, get_interactions2, standardCorr_column, geneName_col, Ids, NetWorkCutThreshold, file= paste("NetWorkCutThreshold_", NetWorkCutThreshold, "_CoexpressionNetworkNetrankfile.RData", sep=""))

  print(paste("Step finished. Netrank can be applied now. The Data has been saved to ", workingDir, "/", datasetDir, "/", paste("NetWorkCutThreshold_", NetWorkCutThreshold, "_CoexpressionNetworkNetrankfile.RData"), sep=""))

}

#' Clusteranalysis on the coexpression network
#' @description This function identifies clusters (modules) in the coexpression network, as well as hub genes.
#' @param datasetDir Name of the directory containing the dataset. The directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param networkType Type of the previously computed coexpression network, default = unsigned.
#' @param power_thre Soft power threshold used for the coexpression network creation.
#' @return Creates a .csv file containing the cluster membership for each gene, a .csv file containing all hub genes, a .csv file containing detailed information regarding the most common GO terms for each cluster, as well as a .csv file containing each genes correlation towards the specific phenotype.
#' @examples
#' check README for an example
#' @export
CoexpressionNetwork_Clusteranalysis <- function(datasetDir, workingDir, networkType = "unsigned", power_thre, numOfGOterms = 10){
  
  setwd(workingDir)
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);

  print("Loading network information")
  loaddata <- load(file = paste(workingDir, datasetDir, '/',networkType,"/CoexpressionNetwork.RData",sep = ""))
  loadWGCNA <- load(file = paste(workingDir, datasetDir, '/CoexpressionNetwork_TOM-block.1.RData', sep = "")) # load WGCNA Network info

  mergedColors = WGCNA::labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  pdf(file = paste(datasetDir, '/',networkType,"/Plots/2_dendrogram.pdf",sep = ""), width = 12, height = 9);
  WGCNA::plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  print("Computing modules (clusters), hub genes and GO term analysis.")
  # module membership of genes
  moduleLabels = net$colors
  moduleColors = WGCNA::labels2colors(net$colors)
  # module eigengenes (1. principal component of modules)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];

  # identify hub genes
  hubs = WGCNA::chooseTopHubInEachModule(datExpr, power = power_thre, colorh=moduleColors, type=networkType)

  nSamples = nrow(datExpr)

  MEs0 = WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes  # calculate the 1st Principal Component (PC) i.e., module eigengene (ME), of each module this time with color labels
  MEs = WGCNA::orderMEs(MEs0)

  #Gene relationship to trait and important modules: Gene Significance and Module Membership
  #We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as
  #(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
  #measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This
  #allows us to quantify the similarity of all genes on the array to every module.


  # somewhat complicated conversion of factor to numeric
  pheno = as.data.frame(as.numeric(levels(datTraits$pheno))[as.integer(datTraits$pheno)]);
  names(pheno) = "pheno"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)

  geneModuleMembership = as.data.frame(WGCNA::cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");

  geneTraitSignificance = as.data.frame(WGCNA::cor(datExpr, pheno, use = "p"));
  GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

  names(geneTraitSignificance) = paste("GS.", names(pheno), sep="");
  names(GSPvalue) = paste("p.GS.", names(pheno), sep="");


  # Annotation
  PROBES = names(datExpr)
  PROBES<- as.character(PROBES)


  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host = 'ensembl.org')
  PROBES2annot <- biomaRt::getBM(attributes = c("affy_hg_u133_plus_2","entrezgene_id", "hgnc_symbol","ensembl_gene_id_version"), mart = mart,   useCache = FALSE,  values=PROBES, filters    = "ensembl_gene_id_version")



  gts=data.frame(geneTraitSignificance, stringsAsFactors=FALSE)
  gts$index=rownames(gts)

  gspv=data.frame(GSPvalue, stringsAsFactors=FALSE)
  gspv$index=rownames(gspv)

  ant=data.frame(PROBES2annot, stringsAsFactors=FALSE)
  ant$index=ant$hgnc_symbol

  genesColors <-data.frame(names(datExpr), stringsAsFactors=FALSE)
  genesColors$colors <-data.frame(moduleColors, stringsAsFactors=FALSE)
  names(genesColors)=c('index', 'color')
  # gene correlation with each modules eigengene
  gmm=data.frame(geneModuleMembership, stringsAsFactors=FALSE)
  gmm$index=rownames(gmm)
  # Student asymptotic p-value for each gene correlation with each modules eigengene
  gmmp=data.frame(MMPvalue, stringsAsFactors=FALSE)
  gmmp$index=rownames(gmmp)

  full_geneInfo <- merge(gts, gspv, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, genesColors, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, gmm, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, gmmp, by = "index", all = TRUE)

  full_geneInfo_annot <- merge(full_geneInfo, ant, by = "index", all = TRUE)
  geneOrder = order(full_geneInfo$color, -abs(full_geneInfo$p.GS.pheno));
  full_geneInfo = full_geneInfo[geneOrder, ]

  full_geneInfo <- merge(gts, gspv, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, genesColors, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, gmm, by = "index", all = TRUE)
  full_geneInfo <- merge(full_geneInfo, gmmp, by = "index", all = TRUE)

  full_geneInfo_annot <- merge(full_geneInfo, ant, by = "index", all = TRUE)

  ## GO Analysis
  #Enrichment analysis directly within R
  #The WGCNA package now contains a function to perform GO enrichment analysis using a simple, single step. To
  #run the function, Biconductor packages GO.db, AnnotationDBI, and the appropriate organism-specific annotation
  #package(s) need to be installed before running this code. The organism-specific packages have names of the form
  #org.Xx.eg.db, where Xx stands for organism code, for example, Mm for mouse, Hs for human, etc. The only exception
  #is yeast, for which no org.Xx.eg.db package is available; instead, the package carries the name org.Sc.sgd.db. Please
  #visit the Bioconductor main page at http://www.bioconductor.org to download and install the required packages.
  #=====================================================================================


  full_geneInfo_annot_unique=full_geneInfo_annot[!duplicated(full_geneInfo_annot$entrezgene_id), ]
  GOenr= WGCNA::GOenrichmentAnalysis(full_geneInfo_annot_unique$color, full_geneInfo_annot_unique$entrezgene_id, organism = "human", nBestP = numOfGOterms);


  #The function runs for awhile and returns a long list, the most interesting component of which is
  tab = GOenr$bestPTerms[[4]]$enrichment

  #This is an enrichment table containing the 10 best terms for each module present in moduleColors. Names of the


  #We refer the reader to the help page of the function within R (available using ?GOenrichmentAnalysis at the R prompt)
  #for details of what each column means. Because the term definitions can be quite long, the table is a bit difficult to
  #display on the screen. For readers who prefer to look at tables in Excel or similar spreadsheet software, it is best to
  #save the table into a file and open it using their favorite tool:

  save(GOenr,full_geneInfo_annot, full_geneInfo_annot_unique, moduleLabels, moduleColors, hubs, file= paste(datasetDir, '/',networkType,"/ClusterAnalysis.RData",sep = ""))
  write.table(hubs, file = paste(datasetDir, '/',networkType,"/HubGenes.csv",sep = ""), row.names = FALSE)
  write.table(data.frame(moduleLabels, moduleColors), file = paste(datasetDir, '/',networkType,"/GeneModules.csv",sep = ""), row.names = T)
  write.table(tab, file = paste(datasetDir, '/',networkType,"/GOEnrichmentTable.csv",sep = ""), quote = T, row.names = F)
  write.table(full_geneInfo_annot_unique[,1:3], file = paste(datasetDir, '/',networkType,"/PhenoCorr.csv",sep = ""), row.names = F)
  print("Cluster and GO Term analysis finished.")
  print(paste("The 10 best GO terms for each module have been saved to ", workingDir, "/", datasetDir, "/", networkType, "/GOEnrichmentTable.csv", sep= ""))
  print(paste("The cluster information for each gene has been saved to ", workingDir, "/", datasetDir, "/", networkType, "/GeneModules.csv", sep= ""))
  print(paste("The hub genes for each module have been saved to ", workingDir, "/", datasetDir, "/", networkType, "/HubGenes.csv", sep= ""))
  print(paste("The phenotype correlation for each gene has been saved to ", workingDir, "/", datasetDir, "/", networkType, "/PhenoCorr.csv", sep= ""))
  print(paste("Objects for further Analsysis have been saved to ", workingDir, "/", datasetDir, "/", networkType, "/ClusterAnalysis.RData", sep= ""))
}

#' Fetches cluster information and hub genes from a gene signature
#' @description This function takes a gene signature and outputs cluster information and hub genes for this signature.
#' @param datasetDir Name of the directory containing the dataset. The directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param networkType Type of the previously computed coexpression network, default = unsigned.
#' @param NetWorkCutThreshold Previously chosen threshold for cutting edges with low coexpression values.
#' @param p_threshold p-value threshold for the gene signature. Genes with a p-value greater or equal to p_treshold are discarded.
#' @param numberOfGenes Number of genes that will be in the signature.
#' @param pathToSignatureFile Path to file that contains a ranked list of genes, i.e. the NetRank output file.
#' @return This function saves three files in a .csv format. The first file contains a signature of genes with full cluster information, the second file lists all hub genes found in the signature, the third file contains additional information about the clusters found in the signature.
#' @examples
#' check the README for an example
#' @export
SignatureClusters <- function(p_threshold, numberOfGenes, pathToSignatureFile, workingDir, datasetDir, rankCol = "rank_NetRank_score", networkType){
  # This program takes the top n Ranked genes with a p-value lower than the threshold and fetches the cluster memberships for those genes.
  # load signature file
  options(stringsAsFactors = FALSE);
  setwd(workingDir)
  GeneSignature <- read.csv(pathToSignatureFile, sep=",")
  hubs <- read.csv(paste(workingDir, datasetDir, "/",networkType, "/HubGenes.csv", sep=""), sep = " ")
  # load cluster info
  ClusterInfo <- load(file = paste(workingDir, datasetDir, "/",networkType, "/ClusterAnalysis.RData", sep=""))
  # merge both into one frame
  GeneSignature <- merge(GeneSignature, full_geneInfo_annot_unique, by.x = "GeneName", by.y = "index")
  # filter by phenotype p value
  GeneSignature = subset(GeneSignature, GeneSignature$p.GS.pheno < p_threshold)
  # oder by NetRank ranking
  GeneSignature = GeneSignature[order(GeneSignature[,rankCol]),]
  # take top n genes
  top_n_genes = head(GeneSignature[order(GeneSignature[,rankCol]),], n=numberOfGenes)
  top_n_genes = merge(data.frame(moduleLabels, moduleColors), top_n_genes, by.x = 0, by.y= "GeneName")
  colnames(top_n_genes)[colnames(top_n_genes) == "Row.names"] <- "GeneName"
  removecolumns <- c("color", "moduleLabels")
  top_n_genes = top_n_genes[, !(names(top_n_genes) %in% removecolumns)]
  top_n_genes = top_n_genes[order(top_n_genes[, rankCol]),]
  
  tab = GOenr$bestPTerms[[4]]$enrichment


  modules = c()
  totalmodulesize = c()
  singatureGenes_inModule = c()
  # get all modules present in gene signature
  for(i in unique(top_n_genes$moduleColors)){
    modules = c(modules, i)
  }
  # get total module size
  for (i in modules){
    # take modSize from tab
    totalmodulesize = c(totalmodulesize, tab[tab$module == i, 2][1])  }
  # get number of signature genes in each module
  for (i in modules){
    singatureGenes_inModule = c(singatureGenes_inModule, nrow(subset(top_n_genes, top_n_genes$moduleColors == i)))
  }
  signatureclusters <- data.frame(modules, totalmodulesize, singatureGenes_inModule)

  hubs = merge(hubs, top_n_genes, by.x = "x", by.y = "GeneName")

  write.table(top_n_genes, file = paste(datasetDir, '/',networkType,"/signature_with_clusters.csv",sep = ""), row.names = FALSE)
  write.table(hubs, file = paste(datasetDir, '/',networkType,"/signature_hubs.csv",sep = ""), row.names = FALSE)
  write.table(signatureclusters, file = paste(datasetDir, '/',networkType,"/signatureclusters.csv",sep = ""), row.names = FALSE)
  print("The clusters and hub genes were saved in the following location:")
  print(paste("The top ",numberOfGenes," from the signature with the full cluster and GO Term information have been saved to:"), sep="")
  print(paste(workingDir, "/", datasetDir, "/", networkType, "/signature_with_clusters.csv", sep= ""))
  print(paste("The hubs found in the signature have been saved to: ", workingDir, "/", datasetDir, "/", networkType, "/signature_hubs.csv", sep= ""))
  print(paste("The clusters found in the signature with total cluster size and  number of genes from the signature in each cluster have been saved to:"), sep="")
  print(paste(workingDir, "/", datasetDir, "/", networkType, "/signatureclusters.csv", sep= ""))



}

