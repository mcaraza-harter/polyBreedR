#' Read SNP array intensity data
#' 
#' Read SNP array intensity data
#' 
#' The first two columns of the tab-delimited input file should be the SNP and Sample ID. Columns labeled "X" and "Y" contain the signal intensities for the two alleles. Use \code{output} to specify whether to return the ratio = Y/(X+Y) or theta = atan(Y/X)*2/pi. 
#' 
#' @param filename filename 
#' @param skip number of lines to skip before the header line with the column names
#' @param output Either "ratio" or "theta"
#' 
#' @return matrix with dimensions markers x individuals
#' 
#' @export
#' @importFrom utils read.table
#' @importFrom tidyr pivot_wider
#' 

readXY <- function(filename,skip,output="ratio") {
  stopifnot(output %in% c("theta","ratio"))
  
  data <- read.table(file=filename, skip=skip, sep="\t", header=T,as.is=T,check.names=F) 
  
  data.1 <- data[,c(colnames(data)[1:2],"X")]
  X <- data.1 %>% pivot_wider(names_from=colnames(data.1)[2],values_from="X")
  #X <- pivot_wider(data=data,id_cols=colnames(data)[1:2],names_from=colnames(data)[2],values_from="X")
  sample.id <- colnames(X)[-1]
  snp.id <- X[,1][[1]]
  X <- as.matrix(X[,-1])
  
  data.2 <- data[,c(colnames(data)[1:2],"Y")]
  Y <- data.2 %>% pivot_wider(names_from=colnames(data.2)[2],values_from="Y")
  #Y <- pivot_wider(data=data,id_cols=colnames(data)[1:2],names_from=colnames(data)[2],values_from="Y")
  Y <- as.matrix(Y[,-1])
  
  if (output=="ratio") {
    data <- Y/(X+Y)
  } else {
    data <- atan(Y/X)*2/3.14159
  }
  rownames(data) <- snp.id
  colnames(data) <- sample.id
  return(data)
}
