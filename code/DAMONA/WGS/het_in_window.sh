#!/bin/bash
#SBATCH --job-name=hetwd
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=batch
#SBATCH --mail-user=nforneris@uliege.be
#SBATCH --mail-type=ALL
#SBATCH --array=1-29%10
#SBATCH --output=hetwd.out
 
 
 # each row in the output file is: snp individual window_heterozygosity 
 
 module load releases/2022b
 module load intel/2022b

 # for 50kb windows
 ifort -O3 -mcmodel=medium het_in_window_50.f90 -o het_in_wd

 # for 25kb windows
 ifort -O3 -mcmodel=medium het_in_window_25.f90 -o het_in_wd

 chromo=$SLURM_ARRAY_TASK_ID

 SCRATCH=$LOCALSCRATCH/$SLURM_JOB_ID
 mkdir -p $SCRATCH
 
 cp '/scratch/ulg/genan/forneris/eval/seq/seq_evalset_gen.txt.gz'      $SCRATCH/.
 cp '/scratch/ulg/genan/forneris/eval/seq/local/het_in_wd/het_in_wd'   $SCRATCH/.
 cd $SCRATCH/

 zcat seq_evalset_gen.txt.gz  | awk -v chr=$chromo '$1==chr' > genfile.txt
 nsnp=`wc -l genfile.txt | awk '{print $1}'`
 echo $nsnp

 ./het_in_wd

 gzip -9 out_het_snp*
 gzip -9 out_het_counts*

 cp out_het_snp_info.txt.gz /scratch/ulg/genan/forneris/eval/seq/local/het_in_wd/out_het_snp_info_$chromo.txt.gz
 cp out_het_counts.txt.gz /scratch/ulg/genan/forneris/eval/seq/local/het_in_wd/out_het_counts_$chromo.txt.gz

 # remove temporary files
 rm -rf $SCRATCH
