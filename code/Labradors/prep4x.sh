#!/bin/bash
#SBATCH --job-name=prep4x
#SBATCH --output=prep4x.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --time=02:00:00

module load releases/2020a
module load BCFtools/1.10.2-GCC-9.3.0
module load VCFtools/0.1.16-GCC-9.3.0

mkdir -p $LOCALSCRATCH/$SLURM_JOB_ID
cd $LOCALSCRATCH/$SLURM_JOB_ID

cp /scratch/ulg/genan/druet/LAB/4X/merged/chr*.vcf.gz* .
bcftools concat c*vcf.gz -Oz -o Labradors4x.vcf.gz

######### perform SNP selection and filtering ###########

vcftools --gzvcf Labradors4x.vcf.gz --out Labradors4x_SNPs --recode --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.05 --recode-INFO-all

bgzip -c Labradors4x_SNPs.recode.vcf > Labradors4x_SNPs.recode.vcf.gz
tabix -p vcf Labradors4x_SNPs.recode.vcf.gz

cp $GLOBALSCRATCH/Labradors/selected_positions.txt .

vcftools --vcf Labradors4x_SNPs.recode.vcf --positions selected_positions.txt --recode --out Labradors4x_SNPs3
bcftools query -f '%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%PL]\n' Labradors4x_SNPs3.recode.vcf > Labradors4x_SNPs4.vcf

######### Reformat in PL format for ZooRoH #####################

cut -f1-5 Labradors4x_SNPs4.vcf > fam.info
cut -f6- Labradors4x_SNPs4.vcf > PL.info

sed -e "s/\./0,0,0/g" PL.info > tmp
sed -e "s/,/\t/g" tmp > tmp2
paste -d "\t" fam.info tmp2 > Labradors4x_SNPs4_PL.txt

######### Save selected positions and copy back the files #######

cp Labradors4x_SNPs4_PL.txt $GLOBALSCRATCH/Labradors/.

cd ..

rm -Rf $SLURM_JOB_ID


