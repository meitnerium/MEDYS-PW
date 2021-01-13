#!/bin/bash
#SBATCH -t 0-48:00
#SBATCH -c 32
#SBATCH --mem=248000

module load imkl/11.1.4.214

export OMP_NUM_THREADS=32
#time /home/dion/workspace/medys-2016/medys > LOG
time /home/dion/medys-2016/medys  > LOG





