#!/bin/bash
#SBATCH --job-name=wgs
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=batch
#SBATCH --array=1-11
#SBATCH --output=wgs.out

module load R

echo $SLURM_JOB_ID

wdir='/scratch/ulg/genan/forneris/eval'

fi="seq_evalset_gen.txt.gz"

SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
mkdir -p $SCRATCH
cd $SCRATCH

cp $wdir/seq/run_zooroh_wgs.R           $SCRATCH/.
cp $wdir/seq/$fi                        $SCRATCH/.

/usr/bin/time -v Rscript run_zooroh_wgs.R $fi $SLURM_ARRAY_TASK_ID

cp *.RData $wdir/seq/.

# remove temporary files
rm -rf $SCRATCH
