#!/bin/bash
#SBATCH --job-name=rrslc
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=batch
#SBATCH --output=rrslc.out

module load R

SCRATCH=$LOCALSCRATCH'/'$SLURM_JOB_ID
mkdir -p $SCRATCH

cp  genome.PstI.1000nt.digest.bed    $SCRATCH/.
cp  RRS_PstI_snp_all.bed             $SCRATCH/.
cp  SNP_seq.bed                      $SCRATCH/.
cp  RRS_panels_generation_lowcov.R   $SCRATCH/.
cp  beagle_chr*.vcf.gz               $SCRATCH/.

cd $SCRATCH

Rscript RRS_panels_generation_lowcov.R

cp *_positions_to_downsample.txt /scratch/ulg/genan/forneris/eval/RRS_lowcov/.
cp reads_to_correct.txt          /scratch/ulg/genan/forneris/eval/RRS_lowcov/.

rm -rf $SCRATCH



