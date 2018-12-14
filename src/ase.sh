
#!/bin/bash

#PBS -l walltime=15:00:00,nodes=1:ppn=2,mem=15gb
#PBS -o /home/hirschc1/lixx5447/projects/biomap/logs
#PBS -e /home/hirschc1/lixx5447/projects/biomap/logs
#PBS -V
#PBS -N ASE181214
#PBS -M lixx5447@umn.edu
#PBS -m abe
#PBS -r n


module load R
Rscript /home/hirschc1/lixx5447/projects/biomap/src/ase.R
