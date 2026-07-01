#!/bin/bash
#SBATCH --job-name=gbs20
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=batch
#SBATCH --output=gbs20.out

module load releases/2022b
module load BCFtools/1.17-GCC-12.2.0

# this filter (-f PASS) is actually not needed as all snps have the same value in that column
bcftools view -f PASS -o filtered_PASS_GBS.vcf.gz fine_tuning_GBS_data-Cow_GBS-250.vcf.gz
bcftools index filtered_PASS_GBS.vcf.gz

bcftools view -v snps -m2 -M2 -r $(seq -s, 1 29) -o filtered_PASS_GBS_biallelic_snps_auto.vcf.gz filtered_PASS_GBS.vcf.gz
bcftools index filtered_PASS_GBS_biallelic_snps_auto.vcf.gz

# gt format
bcftools query -f '%CHROM\t%CHROM"_"%POS\t%POS\t%REF\t%ALT[\t%GT]\n' filtered_PASS_GBS_biallelic_snps_auto.vcf.gz | tr '\t' ' ' | sed 's/"//g' | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' |  sed 's/1\/0/1/g' |  sed 's/1\/1/0/g' | sed 's/\.\/\./9/g'  > filtered_GBS20K_gt.gen

# extract sample ID
bcftools query -l filtered_PASS_GBS_biallelic_snps_auto.vcf.gz > sample_ids_gbs.txt

# sort
sort -k1,1n -k3,3n filtered_GBS20K_gt.gen -o filtered_GBS20K_gt.gen
