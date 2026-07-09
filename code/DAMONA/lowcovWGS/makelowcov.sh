#!/bin/bash
#SBATCH --job-name=for
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=batch
#SBATCH --output=for.out

module load releases/2022b
module load intel/2022b

# modify pcov, seed and files depending on the desired coverage and then compile 
ifort -O3 -mcmodel=medium makelowcov.f90 -o makelowcov

SCRATCH=$LOCALSCRATCH'/'$SLURM_JOB_ID
mkdir -p $SCRATCH

cp Damona132_snp_AD       $SCRATCH/.
cp makelowcov             $SCRATCH/.

cd $SCRATCH
./makelowcov

# modify depending on the desired coverage
cat out_vcf_lcWGS_0.5x | sed 's/\t/ /g' | sed 's/chr//g' | awk '{$1=$1" "$1"_"$2}1' > out_vcf_lcWGS_0.5x_gen.txt

cp  out_* /scratch/ulg/genan/forneris/eval/lowcovWGS/.

rm -rf $SCRATCH
