#!/bin/sh -l 
# FILENAME: submit2.sh 

#SBATCH --job-name=k7_CN_1024
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yu754@purdue.edu
#SBATCH --nodes=1                   
#SBATCH --ntasks=1                   
#SBATCH --account=standby
#SBATCH --cpus-per-task=24
#SBATCH --mem=0
#SBATCH --time=4:00:00
#SBATCH --output=N512_%j.log

# This job's working directory 
module load anaconda
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
echo Running on host `hostname`
echo Time is `date`
echo This jobs runs on the following processors: 

tprocs=1

mpirun -n 8 python3 -u ./main_dynamics.py>& ./logfile