# Starts cluster, export tidyverse
#Author: Greg Ciconetti
# Date: 10/25/2023

goparallel <- function(ncores=7)
{
  cat(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
  cat("\nClosing any open connections...\n")
  closeAllConnections()
  if(exists("cl")) remove(cl)
  cat(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
  cat(paste("\nStarting new cluster with", ncores, "cores...\n"))
  cl <<- parallel::makeCluster(spec = ncores, type="PSOCK")
  cat(" Cluster initiation complete\n")
  cat(paste("\nCurrent Connections: ", dim(showConnections())[1], "\n"))
  cat(paste("\n", exists("cl"), "\n"))
  
  parallel::clusterEvalQ(cl=cl, expr = {require(tidyverse)})
  cat("\n\n***\nThe tidyverse pacakge has been sent to each core.\nDo you need other parallel::clusterEvalQ or parallel::clusterExport calls before running your code?\n****\n")
}