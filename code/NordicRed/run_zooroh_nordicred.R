library(RZooRoH)

sfile        <- 'sample_ids.txt'
genofile     <- "filtered_DAMONA_biallelic_snps_auto_gt.gen"

# other WGS files
#genofile    <- "filtered_PASS_biallelic_snps_auto_gt.gen"
#genofile    <- "filtered_WGS18K_gt.gen"

# other GBS files
#sfile       <- 'sample_ids_gbs.txt'
#genofile    <- "filtered_GBS20K_gt.gen"
#genofile    <- "filtered_GBS18K_gt.gen"

my.data     <- zoodata(genofile, samplefile = sfile, zformat = "gt", min_maf = 0.01)

ccc         <- 1e-2
mycoef      <- rep(ccc,14) 
my.model    <- zoomodel(K = 14, mix_coef = mycoef)
my.resu     <- zoorun(my.model, my.data, localhbd = TRUE, minmix = 0, maxiter = 200)

save.image("zr_nordicred_wgs.RData")

#save.image("zr_nordicred_wgs_lessfiltering.RData")
#save.image("zr_nordicred_WGS18K.RData")
#save.image("zr_nordicred_GBS20K.RData")
#save.image("zr_nordicred_GBS18K.RData")

quit(save="no")
