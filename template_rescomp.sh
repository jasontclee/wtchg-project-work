#!/bin/bash

# Specify a job name (Keep it very short)
#$ -N example-name

# Project name and target queue (Leave as todd.prj*). Long = max 7 days, short = max 1 day. Last letter signifies node (e.g. prjc will enter the C nodes)
#$ -P todd.prjc
#$ -q long.qc

# Run the job in the current working directory
#$ -cwd -j y

# Log locations which are relative to the current
# working directory of the submission
###$ -o output.log
###$ -e error.log

# Parallel environment settings (Number after shmem signals number of cores to be used)
#  For more information on these please see the wiki
#  Allowed settings:
#   shmem
#   mpi
#   node_mpi
#   ramdisk
#$ -pe shmem 4

# Some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "SGE Job ID: $JOB_ID"
echo "SGE Task ID: $SGE_TASK_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Begin writing your script here
# Write a brief description of your code here

cellranger aggr --id=H1 --csv=/well/todd/users/jasonlee/tauseq_data/H1.csv --normalize=mapped
cellranger aggr --id=KO --csv=/well/todd/users/jasonlee/tauseq_data/KO.csv --normalize=mapped

echo $JOB_ID

# End of job script
