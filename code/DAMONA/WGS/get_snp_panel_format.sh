#!/bin/bash
#SBATCH --job-name=fmt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=24G
#SBATCH --partition=batch
#SBATCH --output=fmt.out

module load releases/2022b
module load BCFtools/1.17-GCC-12.2.0

# Extract/create the WGS genotype file from the filtered VCF data (8417679 biallelic SNPs)

# WGS - AD format 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' DAMONA132_CEPE_unrefined.vcf.gz -o Damona132_cepe_comma
cat Damona132_cepe_comma | tr '\t' ' ' | sed 's/,/ /g' > Damona132_snp_AD
rm Damona132_cepe_comma

# WGS - ZooRoH's GT format
zcat DAMONA132_CEPE_refined.vcf.gz \
  | tr '\t' ' ' | sed 's/"//g' \
  | sed 's/0\/0/2/g' \
  | sed 's/0\/1/1/g' \
  | sed 's/1\/0/1/g' \
  | sed 's/1\/1/0/g' \
  | sed 's/0|0/2/g' \
  | sed 's/0|1/1/g' \
  | sed 's/1|0/1/g' \
  | sed 's/1|1/0/g' \
  | sed 's/\.\/\./9/g' \
  | sed 's/\.|\./9/g' \
  | grep -v '^#' | cut -d ' ' -f1,3,2,4,5,10- \
  | awk '{sub(/^chr/, "", $1); print}' \
  | awk '{$1=$1" "$1"_"$2; print}' \
  > seq_evalset_gen.txt
 
# Arrays (medium density - 30K and low density - 6K)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' DAMONA132_CEPE_unrefined.vcf.gz -o Damona132_cepe_gt
cat Damona132_cepe_gt \
  | tr '\t' ' ' | sed 's/"//g' \
  | sed 's/0\/0/2/g' \
  | sed 's/0\/1/1/g' \
  | sed 's/1\/0/1/g' \
  | sed 's/1\/1/0/g' \
  | sed 's/0|0/2/g' \
  | sed 's/0|1/1/g' \
  | sed 's/1|0/1/g' \
  | sed 's/1|1/0/g' \
  | sed 's/\.\/\./9/g' \
  | sed 's/\.|\./9/g' \
  | grep -v '^#' | cut -d ' ' -f1,3,2,4,5,10- \
  | awk '{sub(/^chr/, "", $1); print}' \
  | awk '{$1=$1" "$1"_"$2; print}' \
  > seq_unrefined_gen.txt

# 30K Array - ZooRoH's GT format
awk 'NR==FNR{keep[$1]; next} ($2 in keep)' 30K_ids.txt seq_unrefined_gen.txt > 30K_Array_gen.txt

# 6K Array - ZooRoH's GT format 
awk 'NR==FNR{keep[$1]; next} ($2 in keep)' 6K_ids.txt seq_unrefined_gen.txt > 6K_Array_gen.txt
