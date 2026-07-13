#!/bin/bash
#SBATCH --job-name=seqvcf
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=0-28
#SBATCH --mem-per-cpu=8192
#SBATCH --partition=batch
#SBATCH --output=seqvcf_%a.out
#
module load Java/1.8.0_281

it=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29')

progdir='/CECI/home/ulg/genan/forneris/programs/beagle'
dirhome='/scratch/ulg/genan/forneris/eval/seq/beagle'

SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID

mkdir -p $SCRATCH

i=${it[${SLURM_ARRAY_TASK_ID}]}

seqfile="chr"$i".vcf.gz"

cp $dirhome"/"$seqfile       $SCRATCH/.

cd $SCRATCH
  zcat $seqfile > vcf.vcf

  mkdir -p $SCRATCH'/TMP/'

  java -Xmx7500m -Djava.io.tmpdir=$SCRATCH'/TMP/' -jar $progdir/beagle.25Mar22.4f6.jar seed=050483 nthreads=2 gt=vcf.vcf out=beagle_chr$i
  rm -r TMP

cp $SCRATCH/beagle_chr$i* $dirhome/.

rm -rf $SCRATCH
