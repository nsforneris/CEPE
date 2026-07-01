#!/bin/bash
#SBATCH --job-name=for
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --partition=batch
#SBATCH --output=for.out

module load releases/2022b
module load intel/2022b

# modify pcov, seed and files depending on the desired density/coverage and then compile 
ifort -O3 -mcmodel=medium makeRRSlowcov.f90 -o makeRRSlowcov

SCRATCH=$LOCALSCRATCH'/'$SLURM_JOB_ID
mkdir -p $SCRATCH

cp reads_to_correct.txt                $SCRATCH/.
cp RRS15K_positions_to_downsample.txt  $SCRATCH/. #modify depending on the desired density/coverage
cp Damona132_snp_AD                    $SCRATCH/.

cp makeRRSlowcov                       $SCRATCH/.

cd $SCRATCH
./makeRRSlowcov

# modify depending on the desired density/coverage
cat out_vcf_RRS15K_2x | sed 's/\t/ /g' | sed 's/chr//g' | awk '{$1=$1" "$1"_"$2}1' > out_vcf_RRS15K_2x_gen.txt

cp out_* /scratch/ulg/genan/forneris/eval/rrs_lowcov/.

rm -rf $SCRATCH
