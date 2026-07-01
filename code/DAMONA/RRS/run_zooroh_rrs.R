library(RZooRoH)

genofile <- "RRS15K_gen.txt"
#genofile <- "RRS30K_gen.txt"
#genofile <- "RRS50K_gen.txt"

my.data  <- zoodata(genofile)

ccc = 1e-2
mycoef  = rep(ccc,14)

my.model <- zoomodel(K = 14, mix_coef = mycoef)
my.resu  <- zoorun(my.model, my.data, localhbd = TRUE, nT = 2, minmix = 0)

save.image("zr_RRS15K_Array.RData")
#save.image("zr_RRS30K_Array.RData")
#save.image("zr_RRS50K_Array.RData")

quit(save="no")
