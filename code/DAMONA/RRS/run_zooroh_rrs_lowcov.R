library(RZooRoH)

genofile = "out_vcf_RRS15K_2x_gen.txt"

#genofile = "out_vcf_RRS15K_5x_gen.txt"
#genofile = "out_vcf_RRS30K_2x_gen.txt"
#genofile = "out_vcf_RRS30K_5x_gen.txt"

my.data     <- zoodata(genofile, zformat = "ad", min_maf = 0.01, freqem = TRUE)

ccc = 1e-2
mycoef  = rep(ccc,14)

my.model    <- zoomodel(K = 14, mix_coef = mycoef)
my.resu     <- zoorun(my.model, my.data, localhbd = TRUE, nT = 2, maxiter = 200, minmix = 0)

save.image("zr_RRS15K_2x.RData")

#save.image("zr_RRS15K_5x.RData")
#save.image("zr_RRS30K_2x.RData")
#save.image("zr_RRS30K_5x.RData")

quit(save="no")
