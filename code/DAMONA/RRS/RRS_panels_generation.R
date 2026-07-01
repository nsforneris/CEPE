##############################################################################
# SIMULATION OF RRS DATA FROM WGS DATA
# (RRS15K / RRS30K / RRS50K panels)
##############################################################################
#
# this script:
#
#   1) Takes the in silico digestion of the genome with the PstI enzyme
#      (already done previously with GBSX; here only the resulting
#      fragments are read) and the list of sequence (WGS) SNPs that fall
#      within those fragments (or very close to the enzyme's cut site).
#
#   2) Identifies SNPs located within +/-3 bp of the restriction site.
#      These SNPs can cause "allelic dropout": if an individual is
#      heterozygous at that site, the haplotype carrying the alternative
#      allele is not properly cut/amplified, so in practice only the
#      haplotype with the reference allele is "seen." If the individual is
#      homozygous for the alternative allele at the restriction site, the
#      entire fragment ends up with no genotype (missing) because the RRS
#      fragment is not generated on either haplotype.
#
#   3) Applies this haplotype correction and calculates the resulting
#      "call rate" for the SNPs in those affected fragments.
#
#   4) Filters markers by call rate (0.95) and by MAF (0.01).
#
#   5) From the already-filtered set of RRS SNPs, builds the different
#      panels according to fragment length:
#         - "large" panel (approx. 200-300 bp)          -> ~ RRS50K 
#         - panel of fragments >=250 bp                 -> ~ RRS30K 
#         - 50% subsample of that last panel            -> ~ RRS15K 
#      (the final numbers depend on the input data; the fragment selection
#      criterion is what defines the panel).
#
#   6) Exports the final VCF files
#
##############################################################################

library(dplyr)
library(stringr)

##############################################################################
# 0. INPUT READING: fragments from the in silico digestion and WGS SNPs
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
# WGS SNPs that fall within the PstI fragments (candidate "RRS" SNPs)
snpgbs = read.table('RRS_PstI_snp_all.bed')
snpgbs = snpgbs[,c(1:2)]
colnames(snpgbs)<- c('ch','pos')
snpgbs$name = paste(snpgbs$ch, snpgbs$pos, sep="_")

# All sequence (WGS) SNPs ###
# used to search for variants right at the
# restriction site (outside the fragment but a few bp away from it)
snpseq = read.table('SNP_seq.bed')
snpseq = snpseq[,c(1:2)]
colnames(snpseq)<- c('ch','pos')
snpseq$name = paste(snpseq$ch, snpseq$pos, sep="_")

# The following steps 1 and 2 build a table with the SNPs that are inside
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

# "rem" gathers all the SNPs located at the restriction site (inside or outside the fragment)
# that will trigger the haplotype correction

rem = rbind(remin, remout)
rem = rem[order( rem[,1], rem[,2]),]
rm(f1, f2, f3, fm1, fm2, fm3, remin, remout)

##############################################################################
# 3. Haplotype correction for allelic dropout at the restriction site
##############################################################################
# For each individual, its genotype at
# the restriction-site SNP(s) ("refs") within the fragment is examined:
#
#   - If only ONE haplotype carries the reference allele at that site
#     (heterozygous at the cut site): it is assumed that only that
#     haplotype was amplified, and that haplotype is "copied" onto the
#     other one for the whole fragment (allelic dropout -> homozygosity is
#     forced).
#   - If NEITHER haplotype carries the reference allele (homozygous
#     alternative at the cut site): the fragment is not generated -> the
#     entire fragment is set to missing (".") for that individual.
#   - If BOTH haplotypes carry the reference allele: there is no dropout,
#     the genotype is left as is.
#
# Note: since individuals were previously phased, there are no missing
# input values (that's why ".|." is not searched for in the original VCF).
##############################################################################

gbs = list()
for (ch in 1:29){
  print(ch)
  filename = paste("beagle_chr",ch,".vcf.gz", sep="")
  vcfseq = read.table(filename, skip=9) # skip header
  inds = names(vcfseq)[10:ncol(vcfseq)]
  vcfseq$name = paste(vcfseq$V1,vcfseq$V2,sep="_")
  vcfseq$name = gsub("chr","",vcfseq$name)
  rownames(vcfseq) = vcfseq$name

  # only the candidate RRS SNPs are kept (inside fragments,
  # plus those at -3 bp from the restriction site)`
  
  vcfseq = vcfseq[unique(snpgbs_ex[snpgbs_ex$ch==ch,'name']),]
  vcfseq = vcfseq %>% mutate_all(as.character)
  rownames(vcfseq) = vcfseq$name
  vcfseq2 = vcfseq

  # restriction-site SNPs that determine which haplotype "survives"
  remch = rem[rem$ch==ch,]
  remfragid = unique(remch$frag)
  gc()

  for (f in 1:length(remfragid)){
    # sub-vcf with only the SNP in the fragment
    snpinfrag = snpgbs_ex[snpgbs_ex$frag == remfragid[f],'name']
    vcffrag = vcfseq[ snpinfrag ,]
    # SNPs that mark the reference allele at the restriction site
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
        # only h1 has the ref allele -> h1 is copied over h2
        count0 = count0 + 1
        ind[,c('h2')] = ind[,c('h1')]
        vcffrag2[,ani] = paste(ind$h1,ind$h2,sep="|")
      } else if ( sum(ind[refs,'h1']) != 0 & sum(ind[refs,'h2']) == 0 ) {
        # only h2 has the ref allele -> h2 is copied over h1
        count0 = count0 + 1
        ind[,c('h1')] = ind[,c('h2')]
        vcffrag2[,ani] = paste(ind$h1,ind$h2,sep="|")
      } else if ( sum(ind[refs,'h1']) != 0 & sum(ind[refs,'h2']) != 0  ) {
        # neither haplotype carries the ref allele -> missing fragment
        count1 = count1 + 1
        ind[,c('h1','h2')] = c('.','.')
        vcffrag2[,ani] = paste(ind$h1,ind$h2,sep="|")
      } else {
        # both haplotypes carry the ref allele -> unchanged
        count2 = count2 + 1
        vcffrag2[,ani] = paste(ind$h1,ind$h2,sep="|")
      }
    }
    print (paste(count0," inds with 1 ref hap for frag ", remfragid[f]))
    print (paste(count1," inds with 0 ref hap for frag ", remfragid[f]))
    print (paste(count2," inds with 2 ref hap for frag ", remfragid[f]))

    # call rate of SNP in the fragment, after the correction
    for (k in rownames(vcffrag2)){
      vcffrag2[k,'call'] = sum(vcffrag2[k,inds]!=".|.")/length(inds)
    }

    # the corrected fragment is written back into the chromosome's vcf
    vcfseq2[snpinfrag,c(inds,'call')] = vcffrag2[snpinfrag,c(inds,'call')]
    rm(vcffrag, vcffrag2, snpinfrag, refs)
  }
  # SNPs located right at the restriction site are discarded
  # (they already served their purpose in defining the dropout,
  # and now they are monomorphic / not genotyped as RRS markers)
  vcf = vcfseq2[!(rownames(vcfseq2) %in% remch$name),]
  rm(vcfseq2, vcfseq, remch, remfragid)
  gbs[[ch]] = vcf
  gc()
}

##############################################################################
# PLOT: call rate distribution in affected fragments
##############################################################################

call = gbs[[1]][!is.na(gbs[[1]]$call),c('name','call')]
for (ch in 2:29){
  call = rbind(call,
               gbs[[ch]][!is.na(gbs[[ch]]$call),c('name','call')])
}

tiff("callrate.tif", compression = "lzw")
hist(call$call, breaks=20, xlab = "Call rate", ylab = "Counts", main ="SNP-Call rate for fragments with a SNP in the restriction site")
dev.off()

rm(ch,f,i,k,count0,count1,count2,ani,ind,vcf,filename,call)

##############################################################################
# 4. AF computation
##############################################################################

for (ch in 1:29){
  nsnp = nrow(gbs[[ch]])
  gbs[[ch]][,'freq']=NA
  print(paste(nsnp," snp in fragments of chr ", ch, sep=""))
  for(s in 1:nsnp){
    nam = sum(gbs[[ch]][s,inds]!=".|.")
    if(nam!=0){
      c2 = sum(gbs[[ch]][s,inds]=="0|0")
      c1.0 = sum(gbs[[ch]][s,inds]=="0|1")
      c1.1 = sum(gbs[[ch]][s,inds]=="1|0")
      gbs[[ch]][s,'freq'] = (c2*2 + c1.0 + c1.1)/(2*nam)
    }
  }
}
rm(c1.0,c1.1,c2,s,snpseq)

##############################################################################
# 5. Quality control de calidad: call rate < 0.95 y MAF < 0.01
##############################################################################

gbscall = gbs
# remove SNP with call rate < 0.95 
qc_call = 0
for (ch in 1:29){
  logi = !is.na(gbscall[[ch]]$call) & gbscall[[ch]]$call<0.95
  gbscall[[ch]] = gbscall[[ch]][ !logi  ,]
  qc_call = qc_call + nrow(gbscall[[ch]])
}

gbsmaf = gbscall
# remove SNP with maf < 0.01 
qc_call_maf = 0
for (ch in 1:29){
  gbsmaf[[ch]] = gbsmaf[[ch]][gbsmaf[[ch]]$freq<=0.99,]
  gbsmaf[[ch]] = gbsmaf[[ch]][gbsmaf[[ch]]$freq>=0.01,]
  qc_call_maf = qc_call_maf + nrow(gbsmaf[[ch]])
}

rm(gbscall,gbs,snpgbs_ex)

# This only recalculates, for reference/comparison purposes,
# how many SNPs would have remained if, instead of correcting fragments they had simply been discarded
filename = paste("beagle_chr",ch,".vcf.gz", sep="")
vcfseq = read.table(filename, skip=9)

gbsprev = gbsmaf
# remove the rest of snp in which we computed call rate ####
qc_prev = 0
for (ch in 1:29){
  logi = !is.na(gbsprev[[ch]]$call)
  gbsprev[[ch]] = gbsprev[[ch]][ !logi  ,]
  qc_prev = qc_prev + nrow(gbsprev[[ch]])
}

##############################################################################
# 6. RRS paneles definition based on fragment length (RRS50K / RRS30K / RRS15K)
##############################################################################

# Original panel: all fragments lengths 200-300 pb (aprox. RRS50K)

# Fragments of 250 - 300 bp -> aprox. RRS30K
len250 = unique(frag[frag$len>=250,'name'])
snpgbs250 = snpgbs[snpgbs$frag %in% len250,'name']

# Subsampling of 50% of these fragments -> aprox. RRS15K
len250sub = sample(len250, length(len250)/2)
snpgbs250sub = snpgbs[snpgbs$frag %in% len250sub,'name']

##############################################################################
# 7. Export final VCF
##############################################################################

# use header + fieldnames from original VCF file
system("zcat beagle_chr1.vcf.gz | head -n 8 > first")
system("zcat beagle_chr1.vcf.gz | awk 'NR==9' | awk '{$1=$1};1' > second")

system("cp second ff")
system("cp second ff250")
system("cp second ff250sub")
expo = gbsmaf
for (ch in 1:29){
  expo[[ch]]$ch = gsub("chr","",expo[[ch]]$V1)
  expo[[ch]] = expo[[ch]][  order( as.numeric(expo[[ch]][,'ch']),
                                   as.numeric(expo[[ch]][,'V2']) ),]
  for(line in 1:nrow(expo[[ch]])){
    expo[[ch]][line,] = gsub("|","/",expo[[ch]][line,],fixed=T)
  }
  write.table(file = 'ff', expo[[ch]][,1:(length(inds)+9)],
              row.names=F, col.names=F, quote=F, append=T)
  write.table(file = 'ff250',expo[[ch]][ (expo[[ch]]$name  %in% snpgbs250) ,1:(length(inds)+9)],
              row.names=F, col.names=F, quote=F, append=T)
  write.table(file = 'ff250sub',expo[[ch]][ (expo[[ch]]$name %in% snpgbs250sub) ,1:(length(inds)+9)],
              row.names=F, col.names=F, quote=F, append=T)
}
system("cat ff | sed 's/ /\t/g' > last")
system("cat first last > RRS_PstI_200_300bp.vcf")     #  ~ RRS50K
system("cat ff250 | sed 's/ /\t/g' > last")
system("cat first last > RRS_PstI_250_300bp.vcf")     #  ~ RRS30K
system("cat ff250sub | sed 's/ /\t/g' > last")
system("cat first last > RRS_PstI_250_300bp_sub.vcf") #  ~ RRS15K
