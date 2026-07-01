# An example of how to split RZooRoH results into several runs (convenient when localhbd = TRUE):
# Rscript run_zooroh_lcwgs.R out_vcf_lcWGS_0.5x_gen.txt 0.5 2 

library(RZooRoH)

args      <- commandArgs(trailingOnly=TRUE)
genofile  <- args[1]
cov       <- args[2]
it        <- as.integer(args[3]) # the 2nd subset of samples (13-24)

ini = 12 * (it-1) + 1
fin = ini + 12 - 1
seqid = seq(ini,fin,1)

my.data <- zoodata(genofile, zformat = "ad", min_maf = 0.01, freqem = TRUE)

ccc = 1e-2
mycoef  = rep(ccc,14) 

my.model    <- zoomodel(K = 14, mix_coef = mycoef)
my.resu     <- zoorun(my.model, my.data, localhbd = TRUE, nT = 1, minmix = 0, maxiter = 200, ids = seqid)

save.image(paste0("zr_lcwgs_",cov,"x_",ini,"-",fin,".RData"))

quit(save="no")
