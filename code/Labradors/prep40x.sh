#!/bin/bash
#SBATCH --job-name=prep40x
#SBATCH --output=prep40x.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=02:00:00

module load releases/2020a
module load BCFtools/1.10.2-GCC-9.3.0
module load VCFtools/0.1.16-GCC-9.3.0

mkdir -p $LOCALSCRATCH/$SLURM_JOB_ID
cd $LOCALSCRATCH/$SLURM_JOB_ID

cp $GLOBALSCRATCH/LAB/50X/merged/chr*.vcf.gz* .
bcftools concat c*vcf.gz -Oz -o Labradors50x.vcf.gz

vcftools --gzvcf Labradors50x.vcf.gz --out Labradors50x_SNPs --recode --remove-indels --min-alleles 2 --max-alleles 2 --recode-INFO-all
vcftools --vcf Labradors50x_SNPs.recode.vcf --site-depth --out stats_filtered

bgzip -c Labradors50x_SNPs.recode.vcf > Labradors50x_SNPs.recode.vcf.gz
tabix -p vcf Labradors50x_SNPs.recode.vcf.gz

########## repeat extractation after setting depth range #####

vcftools --gzvcf Labradors50x.vcf.gz --out Labradors50x_SNPs2 --recode --remove-indels --min-alleles 2 --max-alleles 2 --recode-INFO-all --min-meanDP 33 --max-meanDP 48
vcftools --vcf Labradors50x_SNPs2.recode.vcf --site-depth --out stats_filtered2

bgzip -c Labradors50x_SNPs.recode2.vcf > Labradors50x_SNPs2.recode.vcf.gz
tabix -p vcf Labradors50x_SNPs2.recode.vcf.gz

########## select biallelic SNPs from original reference files ######

cp $GLOBALSCRATCH/LAB/biallelic_positions1.txt .

vcftools --vcf Labradors50x_SNPs2.recode.vcf --positions biallelic_positions1.txt --recode --out Labradors50x_SNPs3
bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%PL]\n' Labradors50x_SNPs3.recode.vcf > Labradors50x_SNPs4.vcf

######### Reformat in PL format for ZooRoH #####################

sed -e "s/,/\t/g" Labradors50x_SNPs4.vcf > Labradors50x_SNPs4_PL.txt
sed -e "s/\t/ /g" Labradors50x_SNPs4_PL.txt > tmp

cut -f1-5 Labradors50x_SNPs4.vcf > fam.info
cut -f6- Labradors50x_SNPs4.vcf > PL.info

sed -e "s/\./0,0,0/g" PL.info > tmp
sed -e "s/,/\t/g" tmp > tmp2
paste -d "\t" fam.info tmp2 > Labradors50x_SNPs4_PL.txt

######### Save selected positions and copy back the files #######

cut -f 1,3 fam.info > selected_positions.txt

cp Labradors50x_SNPs4_PL.txt $GLOBALSCRATCH/Labradors/.
cp selected_positions.txt $GLOBALSCRATCH/Labradors/.

cd ..

rm -Rf $SLURM_JOB_ID


