#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=2,mem=40gb
#PBS -o /home/hirschc1/lixx5447/projects/biomap/logs
#PBS -e /home/hirschc1/lixx5447/projects/biomap/logs
#PBS -V
#PBS -N expgeneNumb
#PBS -M lixx5447@umn.edu
#PBS -m abe
#PBS -r n


module load R
Rscript /home/hirschc1/lixx5447/projects/biomap/src/expGeneNumb.R
~

