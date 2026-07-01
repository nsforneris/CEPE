##############################################################################
# PREPARATION OF INPUT FILES TO COMBINE RRS WITH LOW_COVERAGE SEQUENCING
# (RRS15K@2x, RRS15K@5x, RRS30K@2x AND RRS30K@5x)
##############################################################################
#
# this script is the "adapted" version of script RRS_panels_generation.R:
# instead of exporting a final, already-corrected VCF, it generates two text
# files that are later used by the Fortran program makeRRSlowcov.f90 to:

#   a) Know at which which positions (SNP) the read downsampling must be
#      performed (file "RRS[15/30]K_positions_to_downsample.txt").
#
#   b) Know, for each individual and each SNP affected by allelic dropout
#      at the restriction site, how its read counts (AD) must be corrected
#      before performing the downsampling (file "reads_to_correct.txt").
#      This replaces the haplotype correction "0|1 -> 0|0" that was done in
#      script RRS_panels_generation.R, now expressed instead as a correction
#      instruction on the reads (AD) rather than on the phased VCF.
#
# The logic for detecting SNP at the restriction site (sections 1 and 2)
# is identical to that in script RRS_panels_generation.R.
# Only step 3 onward changes: instead of rewriting the corrected genotype
# inside the VCF, the correction instruction is stored (which individual,
# which SNP, and whether the read count of the alternative allele
# [correction_val=0] or of both alleles [correction_val=1] must be set to 0)
# to be applied later to the AD data.
#
##############################################################################

library(dplyr)
library(stringr)

##############################################################################
# 0. INPUT READING: fragments from the in silico digestion and WGS SNP
##############################################################################
# fragments ####
# Fragments resulting from the in silico digestion of the genome with PstI
# (previously performed with GBSX). Each row = one fragment (chromosome,
# start position, end position)

frag = read.table('genome.PstI.1000nt.digest.bed')
colnames(frag) <- c('ch','posini','posend')
nfrag = nrow(frag)
for(i in 1:nfrag){
  frag[i,'chr'] <- gsub("chr_", "", frag[i,'ch'])
}
frag$ch = frag$chr; frag$chr=NULL
frag$len = frag$posend-frag$posini
frag$name = seq(1,nrow(frag))

# RRS SNP ####
# WGS SNP that fall within the PstI fragments (candidate "RRS" SNP)
snpgbs = read.table('RRS_PstI_snp_all.bed')
snpgbs = snpgbs[,c(1:2)]
colnames(snpgbs)<- c('ch','pos')
snpgbs$name = paste(snpgbs$ch, snpgbs$pos, sep="_")

# All WGS SNP ###
# used to search for variants right at the
# restriction site (outside the fragment but a few bp away from it)
snpseq = read.table('SNP_seq.bed')
snpseq = snpseq[,c(1:2)]
colnames(snpseq)<- c('ch','pos')
snpseq$name = paste(snpseq$ch, snpseq$pos, sep="_")

# The following steps 1 and 2 build a table with the SNP that are inside
# a fragment and that may eventually need to be corrected for some
# individuals (because they fall near the enzyme's cut site)

##############################################################################
# 1. SNP inside de fragment, at <=3 pb from one of its extremes
##############################################################################

for (i in 1:nfrag){
  ubi = (snpgbs$ch == frag[i,'ch']) & (snpgbs$pos>=frag[i,'posini']) & (snpgbs$pos<frag[i,'posend'])
  snpgbs[ubi,'frag'] = frag[i,'name']
  snpgbs[ubi,'fragini'] = frag[i,'posini']
  snpgbs[ubi,'fragend'] = frag[i,'posend']
}
snpgbs$distini = snpgbs$pos-snpgbs$fragini
snpgbs$distend = snpgbs$fragend-snpgbs$pos

remin = rbind(snpgbs[snpgbs$distini<3,c('ch','pos','frag','name')],
              snpgbs[snpgbs$distend<=3,c('ch','pos','frag','name')])
rm(ubi)

##############################################################################
# 2. SNP outside the fragment, within a -3 bp window (cut site)
##############################################################################

f1 = frag[,c('ch','posini','name')]
f1$posini = f1$posini - 1
f2 = frag[,c('ch','posini','name')]
f2$posini = f2$posini - 2
f3 = frag[,c('ch','posini','name')]
f3$posini = f3$posini - 3
fm1 = frag[,c('ch','posend','name')]
fm2 = frag[,c('ch','posend','name')]
fm2$posend = fm2$posend + 1
fm3 = frag[,c('ch','posend','name')]
fm3$posend = fm3$posend + 2
colnames(fm1) = colnames(f1)
colnames(fm2) = colnames(f1)
colnames(fm3) = colnames(f1)
remout = rbind(f1, f2, f3, fm1, fm2, fm3)
colnames(remout) = c('ch','pos','frag')
remout[,'name'] = paste(remout[,'ch'], remout[,'pos'], sep="_")
remout = remout[remout$name %in% snpseq$name,]

# such SNP may end up associated with more than one fragment
# (for example, if it is located right between two cut sites):
# all fragment-SNP associations are then kept

snpgbs_ex = rbind(snpgbs[,c('ch','pos','frag','name')], remout)
snpgbs_ex = snpgbs_ex[order(snpgbs_ex[,1], snpgbs_ex[,2]),]

# "rem" gathers all the SNP located at the restriction site (inside or outside the fragment)
# that will trigger the haplotype correction

rem = rbind(remin, remout)
rem = rem[order( rem[,1], rem[,2]),]
rm(f1, f2, f3, fm1, fm2, fm3, remin, remout)

##############################################################################
# 3. Detection of allelic dropout and recording of the corrections to be
#    applied to the read counts (AD), instead of to the phased VCF
##############################################################################
# Same as in RRS_panels_generation.R: since individuals were phased, there are no 
# missing input values (".|." is not searched for in the original vcf).
#
# For each fragment and individual, instead of directly correcting the
# genotype, an instruction is stored in "gbs_mod":
#   correction_val = 0  -> heterozygous at the cut site (dropout of one
#                           allele): in makeRRSlowcov.f90, the read count for the
#                           alternative allele (ad2) of the affected SNP is
#                           set to 0 for that individual.
#   correction_val = 1  -> homozygous alternative at the cut site (the
#                           fragment is not generated): in makeRRSlowcov.f90, the
#                           counts for both alleles (ad1 and ad2) are set
#                           to 0, i.e., missing for that individual/SNP.
#   (if homozygous reference at the cut site, no correction is needed: no
#    instruction is added)
##############################################################################

gbs = list()
gbs_mod = list()

for (ch in 1:29){
  print(ch)
  filename = paste("beagle_chr",ch,".vcf.gz", sep="")
  vcfseq = read.table(filename, skip=9)
  inds   = names(vcfseq)[10:ncol(vcfseq)]
  vcfseq$name = paste(vcfseq$V1,vcfseq$V2,sep="_")
  vcfseq$name = gsub("chr","",vcfseq$name)
  rownames(vcfseq) = vcfseq$name
  
  # only the candidate RRS SNP are kept (those inside fragments  
  # plus those at -3bp from the restriction site)
  vcfseq = vcfseq[unique(snpgbs_ex[snpgbs_ex$ch==ch,'name']),]
  vcfseq = vcfseq %>% mutate_all(as.character)
  rownames(vcfseq) = vcfseq$name
  vcfseq2 = vcfseq
  
  # restriction-site SNP that determine which haplotype "survives"
  remch = rem[rem$ch==ch,]
  remfragid = unique(remch$frag)
  gc()
  
  for (f in 1:length(remfragid)){
    # sub-vcf with only the SNP of this single fragment
    snpinfrag = snpgbs_ex[snpgbs_ex$frag == remfragid[f],'name']
    vcffrag = vcfseq[ snpinfrag ,]
    # SNP(s) that mark the reference allele at the restriction site
    refs = remch[remch$frag == remfragid[f],'name']
    
    vcffrag2 = vcffrag
    count0 = 0
    count1 = 0
    count2 = 0

    for (ani in inds){
      ind = vcffrag[,c('name', ani)]
      rownames(ind) = ind$name
      ind[,c('h1','bar','h2')] = str_split_fixed(ind[, ani],"",3)
      ind[,'h1'] = as.integer(ind[,c('h1')])
      ind[,'h2'] = as.integer(ind[,c('h2')])
      if ( sum(ind[refs,'h1']) == 0 & sum(ind[refs,'h2']) != 0 ){
        # only h1 has the REF allele -> allelic dropout
        count0 = count0 + 1
        gbs_mod[[length(gbs_mod) + 1]] <- cbind(snpinfrag,ani,0)
        
      } else if ( sum(ind[refs,'h1']) != 0 & sum(ind[refs,'h2']) == 0 ) {
        # only h2 has the REF allele -> allelic dropout
        count0 = count0 + 1
        gbs_mod[[length(gbs_mod) + 1]] <- cbind(snpinfrag,ani,0)
        
      } else if ( sum(ind[refs,'h1']) != 0 & sum(ind[refs,'h2']) != 0  ) {
        # neither haplotype carries the REF allele -> missing fragment
        count1 = count1 + 1
        gbs_mod[[length(gbs_mod) + 1]] <- cbind(snpinfrag,ani,1)
        
      } else {
        # both haplotypes carry the REF allele -> no correction
        count2 = count2 + 1
      }
    }
  }
  
  # SNP that were located right at the restriction site are discarded
  # from the RRS positions set (they already served their purpose in
  # detecting the dropout)
  vcf = vcfseq2[!(rownames(vcfseq2) %in% remch$name),1:3]
  rm(vcfseq2, vcfseq, remch, remfragid)
  gbs[[ch]] = vcf
  gc()
}

gbs = rbindlist(gbs)
rm(ch,f,count0,count1,count2,ani,ind,vcf,filename)

##############################################################################
# 4. Definition of the RRS15K / RRS30K (fragments >=250 pb, and
#    50% subsampling of those fragments) and export of the files
#    used by the fortran program
##############################################################################

len250 = unique(frag[frag$len>=250,'name'])
snpgbs250 = snpgbs[snpgbs$frag %in% len250,'name']

# a. RRS30K: positions to downsample
gbs1 = gbs
gbs1[, V3 := paste0(V1,"_",V2)]
gbs1[, V3 := gsub("chr","",V3)]
gbs1[, V4 := gsub("chr","",V1)]
gbs1 = gbs1[V3 %in% snpgbs250]
gbs1[, V4 := as.numeric(V4)]
gbs1[, V2 := as.numeric(V2)]
setorder(gbs1, V4, V2)
gbs1[, V4 := NULL]

# Correction instructions (gbs_mod), restricted to the RRS30K set of SNP
gbs_mod1 <- lapply(gbs_mod, as.data.frame)
gbs_mod1 <- rbindlist(gbs_mod1)
gbs_mod1 = gbs_mod1[snpinfrag %in% gbs1$V3]
gbs_mod1[, ani := gsub("V","",ani)]
gbs_mod1[, ani1 := as.integer(ani)-10+1]
gbs_mod1[, ani := ani1]
gbs_mod1[, ani1 := NULL]

gbs = gbs1; rm(gbs1)
gbs_mod = gbs_mod1; rm(gbs_mod1)

# Files read by the Fortran program:
#   - RRS30K_positions_to_downsample.txt : positions of the RRS30K panel
#   - reads_to_correct.txt               : AD corrections for dropout
write.table(file = 'RRS30K_positions_to_downsample.txt', gbs, row.names=F, col.names=F, quote=F)
write.table(file = 'reads_to_correct.txt', gbs_mod, row.names=F, col.names=F, quote=F)

# b. RRS15K: subsampling of 50% of the RRS30K fragments
set.seed(123)  # set seed for reproducibility
len250sub = sample(len250, length(len250)/2)
snpgbs250sub = snpgbs[snpgbs$frag %in% len250sub,'name']
write.table(file = 'RRS15K_positions_to_downsample.txt', gbs[V3 %in% snpgbs250sub], row.names=F, col.names=F, quote=F)

# (the generated files are next used by the program makeRRSlowcov.f90
# to downsample the reads and generate a RRS set with low coverage)
