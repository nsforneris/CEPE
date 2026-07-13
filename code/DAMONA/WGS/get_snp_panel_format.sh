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
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' Damona132anim.vcf.gz -o Damona132_snp_AD

# WGS - ZooRoH's GT format
bcftools query -f '%CHROM\t%CHROM"_"%POS\t%POS\t%REF\t%ALT[\t%GT]\n' Damona132anim.vcf.gz | tr '\t' ' ' | sed 's/"//g' | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' |  sed 's/1\/0/1/g' |  sed 's/1\/1/0/g' |  sed 's/0|0/2/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/0/g' | sed 's/\.\/\./9/g' |  sed 's/\.|\./9/g'  > seq_evalset_gen.txt

# 30K Array - GT format
awk 'NR==FNR{keep[$1]; next} ($2 in keep)' 30K_ids.txt seq_evalset_gen.txt > 30K_Array_gen.txt

# 6K Array - GT format 
awk 'NR==FNR{keep[$1]; next} ($2 in keep)' 6K_ids.txt seq_evalset_gen.txt > 6K_Array_gen.txt
