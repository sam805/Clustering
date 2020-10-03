#*******************Preprocessing Data ******************
library("dplyr")
library("readxl")
library("naniar")

#*************************************************************
#**************** Function to check infinite values **********
#*************************************************************

is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}

#***********************************************************************************************
#*************************** Data preprocessing for Clustering *********************************
#***********************************************************************************************

Preprocessing <- function(path, is.full.data) {
  tblGL <- read_xlsx(path)
  dfgl <- data.frame(tblGL)
  if (!is.full.data) {
    dfgl2 <- dfgl[dfgl$AuthorsGL == "Brito" | dfgl$AuthorsGL == "Felizola",]
  } else {
    dfgl2 <- dfgl
  }
  #gl.filter <- na.omit(dfgl2)
  invalid.r <- get.invalid.rows(dfgl2)
  if(length(invalid.r > 0)) {
    processed.df <- unique(dfgl2[-invalid.r,])
    #processed.df <- dfgl2[-invalid.r, ]
  } else {
    processed.df <- dfgl2
  }
  results.df <- list()
  # preprocessed.df is the data that can be used for clustering 
  results.df$preprocessed.df <- processed.df
  results.df$fullData <- dfgl2
  return(results.df)
}
#*********************************************************************************

#*************************************************************
#**************** Function to check infinite values **********
#*************************************************************

is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}

#***********************************************************************
#***************************get.invalid.rows****************************
#***********************************************************************


get.invalid.rows <- function(dataf) {
  l = 1
  index <- list()
  um.cols <- c("Liq1DensityGL","GasDensityGL","Liq1ViscosityGL","GasViscosityGL","GasLiqSurTensionGL","PipeDiaGL","PipeIncliGL","PipeRoughGL","Liq1SuperficialVelGL","GasSuperficialVelGL","PressureDropGL")
  dataf.new <- dataf[,um.cols]
  # for (i in 1:dim(dataf.new)[1]) {
  #   colnum <- which(dataf.new[i,] == (-999))
  #   index[[l]] <-c(colnum = ifelse(length(colnum) > 0, colnum, 0),rownum = l)
  #   l = l + 1
  # } # END OF for
  #index.df <- as.data.frame(do.call(rbind, index))
  #invalid.rows <- index.df$rownum[index.df$colnum != 0]
  res <- apply(dataf.new, MARGIN = 2, function(x) which(x ==-999))
  invalid.rows <- c(unique(unlist(res)),which(dataf.new$PressureDropGL == 0))
  return (invalid.rows)
} 



