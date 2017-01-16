#!/bin/bash -l
#SBATCH --job-name=cuda
#SBATCH --comment="This is a multithreaded, single GPU job template"
#SBATCH --partition=cudaq
#SBATCH --mail-type=FAIL,TIME_LIMIT_80
##SBATCH --mail-user=zhengguangzhao2014@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=100000
#SBATCH --time=10:00:00
#SBATCH --output=out-%j.out
#SBATCH --error=out-%j.err
#SBATCH --workdir=.

echo "We have $SLURM_JOB_CPUS_PER_NODE core(s) on these nodes:"
# (could also check $SLURM_JOB_NODELIST)
echo "pdf_961geo.py is running[Waiting...]"

python pdf_961geo.py
