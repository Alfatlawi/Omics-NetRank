#' Building a protein interaction network.
#' @description Creates a protein interaction network from correlation data based on STRINGdb protein interaction information.
#' @param datasetDir Name of the directory containing the dataset. The directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param standardCorr_column Name of the column containing the correlation information in the correlation file. Default = 'PearsonCorrelation'.
#' @param geneName_col Name of the column containing the gene names in the correlation file. Default = 'GeneName'
#' @param dataFolder Name of the directory containing the correlation File. The results will be saved into this directory. Default = 'standardGeneScreeningResults'.
#' @param species STRINGdb ID of the species which the data is from. Default = 9606 for human.
#' @param corrFilename Name of the file containing the gene correlation information.
#' @return Saves files for NetRank application in the "datasetDir" directory.
#' @examples check README for an example
#' @export
fetchStringDB <- function(datasetDir, workingDir, standardCorr_column='PearsonCorrelation', geneName_col= "GeneName", dataFolder='standardGeneScreeningResults', corrFilename, species=9606, score_threshold=200){


  setwd(paste(workingDir, datasetDir, sep='' ));

  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);

  ### Read the standard_corr from previous results
  print("Loading data from correlation file.")
  standard_corr = read.csv(paste(workingDir, datasetDir, "/", dataFolder, "/", corrFilename, sep = ""), sep=",");
  standard_corr <- standard_corr[complete.cases(standard_corr[, standardCorr_column]),] # remove rows with missing or Na value in standardCorr_column

  standard_corr[standardCorr_column] = abs(standard_corr[standardCorr_column]) # Sort by correlation

  ### import string_db ### save object string_db and standard_corr locally
  print("Importing StringDB data.")
  string_db <- STRINGdb::STRINGdb$new(version="10", species=species, score_threshold=score_threshold, input_directory="") 

  ### Map each gene in the standard_corr with its ID from StringDB
  print("Mapping StringDB data.")
  standard_corr <- string_db$map( standard_corr,geneName_col, removeUnmappedRows = TRUE)

  ### Remove the duplicates genes and keep the one that has high correlation as they are already sorted.
  standard_corr = standard_corr[!duplicated(standard_corr$STRING_id), ] # remove duplicated

  # get_interactions from stringDB
  print("Fetching interactions from StringDB")
  get_interactions <- string_db$get_interactions(c(standard_corr$STRING_id))

  ## Normalize to have the scores below 1
  get_interactions$combined_score_norm <- get_interactions$combined_score / max(get_interactions$combined_score)

  # this part removes double entries from string_db call
  get_interactions = get_interactions[!duplicated(cbind(get_interactions$from,get_interactions$to)),  ] # remove duplicated

  ### duplciates the genes to have all genes in the column from.
  get_interactions2 = get_interactions
  get_interactions2$from = get_interactions$to
  get_interactions2$to = get_interactions$from
  get_interactions2 = rbind(get_interactions, get_interactions2)
  get_interactions2 = subset(get_interactions2, is.element(get_interactions2$from, standard_corr$STRING_id))
  get_interactions2 = subset(get_interactions2, is.element(get_interactions2$to, standard_corr$STRING_id))

  ## Take the unique ids and their counts
  Ids = unique(get_interactions2$to)
  N = length(Ids)

  print("Aggregating connectivity scores")
  ### sum of connectivity (outputs) of each genes.
  outs <- data.frame(matrix(ncol = 0, nrow = N))
  outs$from = Ids
  outsCount = aggregate(combined_score_norm ~ from, get_interactions2, sum) # sum of connectivity


  outs = merge(outs, outsCount, all.x=T)
  outs[is.na(outs[,2]),2] <- 0.000000001
  colnames(outs) = c('STRING_id','outs_conn')

  ## Merge the sum as a column
  standard_corr = merge(standard_corr, outs, by='STRING_id')

  save(standard_corr, get_interactions2, standardCorr_column, Ids, file="StringDBNetwork.RData")

  print(paste("Results were saved in ", workingDir, datasetDir, "/", "stringDBNetwork.Rdata", sep=""))

}

