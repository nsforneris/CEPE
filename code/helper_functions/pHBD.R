# how to compute inbreeding using HBD classes upto a threshold T and the HBD probabilities 
# uses the cumhbd function in the ZooRoH package
# usage: Rscript pHBD.R 'zr_6K_Array.RData'

library(data.table)
library(RZooRoH)

args <- commandArgs(trailingOnly = TRUE)
zoorohfile <- args[1]

load(zoorohfile)

Tlab <- c(2,4,8,16,32,64,
          128,256,512,1024,2048,
          4096,8192,16384)

pHBD<- data.table(id = my.resu@sampleids)
for (i in seq_along(Tlab)) {
  colname <- paste0("T",Tlab[i])
  pHBD[, (colname) := cumhbd(my.resu, T = as.numeric(Tlab[i]))]
}

pHBD[]
	 
save(pHBD, file = 'pHBD.RData')
