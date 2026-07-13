# An example of how to split RZooRoH results into several runs (convenient when localhbd = TRUE):
# Rscript run_zooroh_wgs.R out_vcf_WGS_gen.txt 2 

library(RZooRoH)

args      <- commandArgs(trailingOnly=TRUE)
genofile  <- args[1]
it        <- as.integer(args[2]) # the 2nd subset of samples (13-24)

ini = 12 * (it-1) + 1
fin = ini + 12 - 1
seqid = seq(ini,fin,1)

my.data <- zoodata(genofile)

ccc = 1e-2
mycoef  = rep(ccc,14) 

my.model    <- zoomodel(K = 14, mix_coef = mycoef)
my.resu     <- zoorun(my.model, my.data, localhbd = TRUE, minmix = 0, ids = seqid)

save.image(paste0("zr_wgs_",ini,"-",fin,".RData"))

quit(save="no")
