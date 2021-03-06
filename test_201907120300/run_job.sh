#!/bin/bash
#$ -cwd -V
#$ -l h_rt=8:00:00
#$ -pe smp 1
#$ -l h_vmem=24G

# print the date:
date

# module loadings:
if [ -r /nobackup/cemac/cemac.sh ] ; then
  . /nobackup/cemac/cemac.sh
fi
module purge
module load user flexpart
module list

# run flexpart:
FLEXPART

# print the date:
date
