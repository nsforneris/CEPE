#!/bin/bash
#SBATCH --job-name=rrs
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=batch
#SBATCH --output=rrs.out

module load R

SCRATCH=$LOCALSCRATCH'/'$SLURM_JOB_ID
mkdir -p $SCRATCH

cp  genome.PstI.1000nt.digest.bed $SCRATCH/.
cp  RRS_PstI_snp_all.bed          $SCRATCH/.
cp  SNP_seq.bed                   $SCRATCH/.
cp  RRS_panels_generation.R       $SCRATCH/.
cp  beagle_chr*.vcf.gz            $SCRATCH/.

cd $SCRATCH

Rscript RRS_panels_generation.R

vcffile='RRS_PstI_250_300bp_sub.vcf'
name='RRS15K'
cat $vcffile | grep -v "#" | awk '{$6="";$7="";$8="";$9=""}1' | tr -s " " | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/0/g' | sed 's/\.\/\./9/g' | sed 's/chr//'| awk '{t=$2; $2=$3;$3=t; print;}' > $name'_gen.txt'

vcffile='RRS_PstI_250_300bp.vcf'
name='RRS30K'
cat $vcffile | grep -v "#" | awk '{$6="";$7="";$8="";$9=""}1' | tr -s " " | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/0/g' | sed 's/\.\/\./9/g' | sed 's/chr//'| awk '{t=$2; $2=$3;$3=t; print;}' > $name'_gen.txt'

vcffile='RRS_PstI_200_300bp.vcf'
name='RRS50K'
cat $vcffile | grep -v "#" | awk '{$6="";$7="";$8="";$9=""}1' | tr -s " " | sed 's/0\/0/2/g' | sed 's/0\/1/1/g' | sed 's/1\/0/1/g' | sed 's/1\/1/0/g' | sed 's/\.\/\./9/g' | sed 's/chr//'| awk '{t=$2; $2=$3;$3=t; print;}' > $name'_gen.txt'

cp *_gen.txt /scratch/ulg/genan/forneris/eval/RRS/.

rm -rf $SCRATCH
