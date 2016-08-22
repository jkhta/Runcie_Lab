#! /bin/bash -l
#SBATCH -D /home/jkta/projects/marinara
#SBATCH -o /home/jkta/projects/marinara/logs/output/out-%j.txt
#SBATCH -e /home/jkta/projects/marinara/logs/error/error-%j.txt
#SBATCH -J marinara

module load R

R CMD BATCH --no-save --no-restore "--args ${SLURM_ARRAY_TASK_ID}" marinara_prop_comp_FA_inh_batch.R logs/Rout-${SLURM_ARRAY_TASK_ID}.txt
echo "Done"