#! /bin/bash -l
#SBATCH -D /home/jkta/projects/jaeger
#SBATCH -o /home/jkta/projects/jaeger/logs/output/out-%j.txt
#SBATCH -e /home/jkta/projects/jaeger/logs/error/error-%j.txt
#SBATCH -J jaeger_opt

module load R

R CMD BATCH --no-save --no-restore "--args ${SLURM_ARRAY_TASK_ID}" optimization_batch_high.R logs/Rout-${SLURM_ARRAY_TASK_ID}.txt
echo "Done"