#' NetRank Algorithm
#'
#' @description Applies the NetRank algorithm to a previously constructed network.
#' @param networkType Choose between two network types "StringDB" or "Coexpression".
#' @param datasetDir Name of the directory containing the dataset. The directory containing the dataset has to be in the working directory 'workingDir'.
#' @param workingDir Path to the working directory.
#' @param cores Number of cores for parallelization.
#' @param d Dampening factor between 0 and 1. This factor weights the two parts of the NetRank formula against each other (node correlation vs. node connectivity). A d value close to 0 increases the importance of node correlation, a value close to 1 increases the importance of node connectivity. For further information check https://doi.org/10.3389/fbinf.2022.780229. Default = 0.5
#' @param minError Threshold for stopping the main Netrank computation as soon as the root mean square deviation of the netrank score between iterations is below minError. Default = 0.01
#' @param dataFolder Name of the directory containing the correlation File. The results will be saved into this directory. Default = 'standardGeneScreeningResults'.
#' @param NetWorkCutThreshold Previously chosen threshold for cutting edges with low coexpression values. Only relevant if network == "Coexpression".
#' @return Saves a file containing gene NetRank score and ranking, as well as the phenotype correlation. Results are saved in the 'dataFolder' in csv format.
#' @examples The Working directory has to contain a directory named 'datasetDir'. The 'datasetDir' directory has to contain a directory named 'dataFolder'.
#' @export
RunNetRank <- function(network, datasetDir, workingDir, cores, d = 0.5, minError=0.01, dataFolder='standardGeneScreeningResults', NetWorkCutThreshold){

if (network == "StringDB") {

  options(stringsAsFactors = FALSE);
  start.time.total <- Sys.time()

  loadData <- load(file = paste(workingDir, datasetDir,'/StringDBNetwork.RData',sep = ""))
  rmsd_vector <- data.frame(matrix(ncol = 2, nrow = 0))
  N = length(Ids)

  ## set initial rmsd to 1
  rmsd=1
  count=0

  ## use the correlation as initial network values (NetRank_Score current step values, NetRank_prev previous step values)
  standard_corr["NetRank_prev"] = standard_corr[standardCorr_column]
  standard_corr["NetRank_score"] = standard_corr[standardCorr_column]

  ## as_FBM is used for parallelization, the NetRank computations will be stored in this matrices
  mat = bigstatsr::as_FBM(standard_corr["NetRank_prev"])
  mat_part1 = bigstatsr::as_FBM(standard_corr["NetRank_prev"])
  mat_part2 = bigstatsr::as_FBM(standard_corr["NetRank_prev"])

  ### initialize cluster for parallelization
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  print('Computing NetRank Update. This will take a while.')
  `%dopar%` <- foreach::`%dopar%`

  while (rmsd > minError) {
    # measure runtime
    start.time <- Sys.time()
    # iterate over every node (i)
    NetRankUpdate <- foreach::foreach (ix=1:length(Ids), .combine = 'c') %dopar% {

      i = Ids[ix] # current node
      ## get information about all nodes connected to node i
      tmpset = subset(get_interactions2, (get_interactions2$to==i)) ## get connectivity information of all nodes that input this node.
      from_ids = tmpset$from ## get the id's of nodes that connect to the node i.
      connected_nodes = subset(standard_corr, is.element(standard_corr$STRING_id, from_ids)) # bring correlation of the required node. Subset the standard_corr dataframe to only contain nodes connected to node i
      connected_nodes = merge(connected_nodes, tmpset, by.x='STRING_id', by.y='from')
      netrank_quotient = (connected_nodes$NetRank_score/(connected_nodes$outs_conn)) ## divide the significance score (correlation in the first iteration) by the sum of outputs of that node.

      ## Compute NetRank Node Update according to Formula
      part1 = ((1-d) * standard_corr[standard_corr$STRING_id==i, standardCorr_column])  # left part of formula: (1-d)*s
      part2 =  d * (sum(connected_nodes$combined_score_norm * netrank_quotient)) # the second part of the formula

      ## On very rare occasions part2 is bigger than one, in this case map it to one.
      if (part2 > 1){
        part2 <- 1
      }

      ## store results in matrices
      mat[standard_corr$STRING_id==i,] = part1 + part2

      mat_part1[standard_corr$STRING_id==i,] = part1
      mat_part2[standard_corr$STRING_id==i,] = part2 # vector of all p2 values #plot it (see how many are >1), interquartile range (point to median)?

      NULL
    }

    ## update columns
    standard_corr["NetRank_prev"] = standard_corr["NetRank_score"]
    standard_corr["NetRank_score"] = mat[]

    ## compute new rmsd
    rmsd = standard_corr["NetRank_score"] - standard_corr["NetRank_prev"]
    rmsd = rmsd* rmsd
    rmsd = sum(rmsd) / N
    rmsd = sqrt(rmsd)

    count = count+1
    end.time <- Sys.time()
    time.taken <- difftime(end.time, start.time, units = "mins")
    print(paste0('The ',count,'. iteration of Netrank is complete. The rmsd is ', rmsd,'. Time taken for this iteration: ',time.taken))

    rmsd_vector = c()


    rmsd_vector = rbind(rmsd_vector,c(count,rmsd))


  }
  standard_corr$mat_part1=mat_part1[]
  standard_corr$mat_part2=mat_part2[]

  parallel::stopCluster(cl)
  colnames(rmsd_vector) = c('iteration','rmsd')

  end.time.total <- Sys.time()
  time.taken.total  <- difftime(end.time.total, start.time.total, units = "mins")
  print(paste0('Computation finished! Total time taken: ',time.taken.total))

  ## rank by correlation
  standard_corr$rank_std <- NA
  standard_corr$rank_std[order(-standard_corr[standardCorr_column])] <- 1:nrow(standard_corr)

  ## rank by NetRank result
  standard_corr$rank_NetRank_score <- NA
  standard_corr$rank_NetRank_score[order(-standard_corr['NetRank_score'])] <- 1:nrow(standard_corr)
  standard_corr <- standard_corr[, c('STRING_id', 'GeneName', standardCorr_column, 'NetRank_score', 'rank_NetRank_score')]


  ## store results in csv
  write.csv(standard_corr,paste(workingDir, datasetDir,"/",dataFolder, '/netRankResults_stringDB_dampF_', d,'_dataset_',datasetDir,"_netRank_StringDBNetwork.csv",sep = ""),row.names=FALSE)
  write.csv(rmsd_vector,paste(workingDir, datasetDir,"/",dataFolder, '/netRankProcess_stringDB_dampF_', d,'_dataset_',datasetDir,"_netRank_StringDBNetwork.csv",sep = ""),row.names=FALSE)


  print('Results are saved!')
}
else if (network == "Coexpression"){

  options(stringsAsFactors = FALSE);
  start.time.total <- Sys.time()

  loadData <- load(file = paste(workingDir, datasetDir, "/",paste("NetWorkCutThreshold_", NetWorkCutThreshold, "_CoexpressionNetworkNetrankfile.RData", sep=""),sep = ""))
  rmsd_vector <- data.frame(matrix(ncol = 2, nrow = 0))
  N = length(Ids)

  ## set initial rmsd to 1
  rmsd=1
  count=0

  ## use the correlation as initial network values (NetRank_Score current step values, NetRank_prev previous step values)
  standard_corr["NetRank_prev"] = standard_corr[standardCorr_column]
  standard_corr["NetRank_score"] = standard_corr[standardCorr_column]

  ## as_FBM is used for parallelization, the NetRank computations will be stored in this matrices
  mat = bigstatsr::as_FBM(standard_corr["NetRank_prev"])
  mat_part1 = bigstatsr::as_FBM(standard_corr["NetRank_prev"])
  mat_part2 = bigstatsr::as_FBM(standard_corr["NetRank_prev"])

  ### initialize cluster for parallelization
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  print('Computing NetRank Update. This will take a while.')
  `%dopar%` <- foreach::`%dopar%`

  while (rmsd > minError) {
    # measure runtime
    start.time <- Sys.time()
    # iterate over every node (i)
    NetRankUpdate <- foreach::foreach (ix=1:length(Ids), .combine = 'c') %dopar% {
      ## generate new file here to check if code is here if print wasnt shown ##
      i = Ids[ix] # current node
      ## get information about all nodes connected to node i
      tmpset = subset(get_interactions2,  ((get_interactions2$to==i) | (get_interactions2$from==i))) ## get connectivity information of all nodes that input this node.

      tmpset$Connected_nodeId=0
      for (tmpset_row in 1:nrow(tmpset)) {
        if (tmpset[tmpset_row,'from']!=i) {tmpset[tmpset_row,'Connected_nodeId']=tmpset[tmpset_row,'from']}
        else {tmpset[tmpset_row,'Connected_nodeId']=tmpset[tmpset_row,'to']}}
      from_ids = unique(tmpset$Connected_nodeId) ## get the id's of nodes that connect to the node i.
      tmpset=tmpset[,c('Connected_nodeId','corr')]
      connected_nodes = subset(standard_corr, is.element(standard_corr$GeneName, from_ids)) # bring correlation of the required node. Subset the standard_corr dataframe to only contain nodes connected to node i
      connected_nodes = merge(connected_nodes, tmpset, by.x=geneName_col, by.y='Connected_nodeId')
      netrank_quotient = (connected_nodes$NetRank_score/(connected_nodes$outs_conn)) ## divide the significance score (correlation in the first iteration) by the sum of outputs of that node.

      ## Compute NetRank Node Update according to Formula
      part1 = ((1-d) * standard_corr[standard_corr$GeneName==i, standardCorr_column])  # left part of formula: (1-d)*s
      part2 =  d * (sum(connected_nodes$corr * netrank_quotient)) # the second part of the formula

      ## On very rare occasions part2 is bigger than one, in this case map it to one.
      if (part2 > 1){
        part2 <- 1
      }

      ## store results in matrices
      mat[standard_corr$GeneName==i,] = part1 + part2
      mat_part1[standard_corr$GeneName==i,] = part1
      mat_part2[standard_corr$GeneName==i,] = part2 # vector of all p2 values #plot it (see how many are >1), interquartile range (point to median)?

      NULL
    }

    ## update columns
    standard_corr["NetRank_prev"] = standard_corr["NetRank_score"]
    standard_corr["NetRank_score"] = mat[]

    ## compute new rmsd
    rmsd = standard_corr["NetRank_score"] - standard_corr["NetRank_prev"]
    rmsd = rmsd* rmsd
    rmsd = sum(rmsd) / N
    rmsd = sqrt(rmsd)

    count = count+1
    end.time <- Sys.time()
    time.taken <- difftime(end.time, start.time, units = "mins")
    print(paste0('The ',count,'. iteration of Netrank is complete. The rmsd is ', rmsd,'. Time taken for this iteration: ',time.taken))

    rmsd_vector = c()
    rmsd_vector = rbind(rmsd_vector,c(count,rmsd))


  }
  standard_corr$mat_part1=mat_part1[]
  standard_corr$mat_part2=mat_part2[]

  parallel::stopCluster(cl)
  colnames(rmsd_vector) = c('iteration','rmsd')

  end.time.total <- Sys.time()
  time.taken.total  <- difftime(end.time.total, start.time.total, units = "mins")
  print(paste0('Computation finished! Total time taken: ',time.taken.total))

  ## rank by correlation
  standard_corr$rank_std <- NA
  standard_corr$rank_std[order(-standard_corr[standardCorr_column])] <- 1:nrow(standard_corr)

  ## rank by NetRank result
  standard_corr$rank_NetRank_score <- NA
  standard_corr$rank_NetRank_score[order(-standard_corr['NetRank_score'])] <- 1:nrow(standard_corr)

  standard_corr <- standard_corr[, c('GeneName', standardCorr_column, 'NetRank_score', 'rank_NetRank_score')]

  ## store results in csv
  write.csv(standard_corr,paste(workingDir, datasetDir,"/",dataFolder, '/netRankResults_Coex_dampF_', d,'_dataset_',datasetDir, '_CutThreshold_', NetWorkCutThreshold, "_netRank_CoexprNetwork.csv",sep = ""),row.names=FALSE)
  write.csv(rmsd_vector,paste(workingDir, datasetDir,"/",dataFolder, '/netRankProcess_Coex_dampF_', d,'_dataset_',datasetDir, '_CutThreshold_', NetWorkCutThreshold,"_netRank_CoexprNetwork.csv",sep = ""),row.names=FALSE)


  print('Results are saved!')
}
else {

  print("Invalid Network Type, specify between 'StringDB' or 'Coexpression' in function call.")
}



}
