#!/bin/bash
#SBATCH --job-name=pass
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=batch
#SBATCH --output=pass.out

module load releases/2022b
module load BCFtools/1.17-GCC-12.2.0

bcftools view -f PASS -o filtered_PASS.vcf.gz fine_tuning_GBS_data-Cow_WGS.vcf.gz
bcftools index filtered_PASS.vcf.gz

bcftools view -v snps -m2 -M2 -r $(seq -s, 1 29) -o filtered_PASS_biallelic_snps_auto.vcf.gz filtered_PASS.vcf.gz
bcftools index filtered_PASS_biallelic_snps_auto.vcf.gz

# gt format
bcftools query -f '%CHROM\t%CHROM"_"%POS\t%POS\t%REF\t%ALT[\t%GT]\n' filtered_PASS_biallelic_snps_auto.vcf.gz | tr '\t' ' ' | sed 's/"//g' | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' |  sed 's/1\/0/1/g' |  sed 's/1\/1/0/g' |  sed 's/0|0/2/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/0/g' | sed 's/\.\/\./9/g' |  sed 's/\.|\./9/g'  > filtered_PASS_biallelic_snps_auto_gt.gen

# extract sample ID
bcftools query -l filtered_PASS_biallelic_snps_auto.vcf.gz > sample_ids.txt

# sorting
sort -k1,1n -k3,3n filtered_PASS_biallelic_snps_auto_gt.gen -o filtered_PASS_biallelic_snps_auto_gt.gen
