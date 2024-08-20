#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb
module load anaconda3/personal
echo "R is about to run"
cp $HOME/Project_2024_main.R $TMPDIR
R --vanilla < $HOME/HPC_project_2024.R
mv HPCoutput_* $HOME/output_files/
echo "R has finished running"
