# how to compute inbreeding using HBD classes upto a threshold T using the Viterbi approach 

# ussage: Rscript pvit.R 'zr_6K_Array.RData'

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
zoorohfile <- args[1]

load(zoorohfile)

Tord <- c('2'=1,'4'=2,'8'=3,'16'=4,'32'=5,'64'=6,
          '128'=7,'256'=8,'512'=9,'1024'=10,'2048'=11,
          '4096'=12,'8192'=13,'16384'=14)

# map length
map <- data.table(my.data@chrbound)
map[,len := my.data@bp[V2]-my.data@bp[V1]]

df <- data.table(my.resu@hbdseg)

# all possible HBD classes from our definition
all_classes <- seq_along(Tord)   # 1:14

# genome length
genome_length <- sum(map$len)

# proportion at each HBD class
pclass <- df[, .(p = sum(length) / genome_length),
             by = .(id, HBDclass)]


# add missing IDs and missing classes
grid <- CJ(id = my.resu@ids,
           HBDclass = all_classes)
pclass <- merge(grid, pclass,
                by = c("id", "HBDclass"),
                all.x = TRUE)

# replace missing proportions by zero
pclass[is.na(p), p := 0]

# cumulative sum over classes
setorder(pclass, id, HBDclass)

pvit_long <- pclass[, .(
  pvit = cumsum(p),
  HBDclass = HBDclass
),
by = .(id)]


# reshape to one row per individual
pvit <- dcast(pvit_long,
              id ~ HBDclass,
              value.var = "pvit",
              fill = 0)

# rename columns
setnames(pvit,
         old = names(pvit)[-1],
         new = names(Tord))
		 
save(pvit, file = 'pvit.RData')
