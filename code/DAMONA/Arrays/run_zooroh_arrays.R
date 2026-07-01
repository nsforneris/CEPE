library(RZooRoH)

genofile <- "6K_Array_gen.txt"
#genofile <- "30K_Array_gen.txt"

my.data  <- zoodata(genofile)

ccc = 1e-2
mycoef  = rep(ccc,14)

my.model <- zoomodel(K = 14, mix_coef = mycoef)
my.resu  <- zoorun(my.model, my.data, localhbd = TRUE, nT = 2, minmix = 0)

save.image("zr_6K_Array.RData")
#save.image("zr_30K_Array.RData")

quit(save="no")
