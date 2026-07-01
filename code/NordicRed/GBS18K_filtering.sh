#!/bin/bash
#SBATCH --job-name=gbs18
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=batch
#SBATCH --output=gbs18.out

module load releases/2022b
module load BCFtools/1.17-GCC-12.2.0

 cat filtered_WGS18K_gt.gen | sort -k1,1n -k3,3n | awk '{print $1 "\t" $3}' > fixed_snps_wgs18k.txt

# keep the clean GBS snps that are also present in the WGS 
 bcftools view -R fixed_snps_wgs18k.txt -o filtered_gbs_in_wgs.vcf.gz filtered_PASS_GBS_biallelic_snps_auto.vcf.gz
 bcftools index filtered_gbs_in_wgs.vcf.gz

 bcftools view -v snps -m2 -M2 -r $(seq -s, 1 29) -o filtered_gbs_in_wgs_biallelic_snps_auto.vcf.gz filtered_gbs_in_wgs.vcf.gz
 bcftools index filtered_gbs_in_wgs_biallelic_snps_auto.vcf.gz

# gt format
 bcftools query -f '%CHROM\t%CHROM"_"%POS\t%POS\t%REF\t%ALT[\t%GT]\n' filtered_gbs_in_wgs_biallelic_snps_auto.vcf.gz | tr '\t' ' ' | sed 's/"//g' | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' |  sed 's/1\/0/1/g' |  sed 's/1\/1/0/g' |  sed 's/0|0/2/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/0/g' | sed 's/\.\/\./9/g' | sed 's/\.|\./9/g' > filtered_GBS18K_gt.gen

# extract sample ID
 bcftools query -l filtered_gbs_in_wgs_biallelic_snps_auto.vcf.gz > sample_ids_gbs.txt

# sorting
 sort -k1,1n -k3,3n filtered_GBS18K_gt.gen -o filtered_GBS18K_gt.gen
