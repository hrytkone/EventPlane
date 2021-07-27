#!/bin/bash
## sbatch usage: sbatch submit.sh <comment> <inputdir> <bmin> <bmax>
## sbatch will check arguments from the comments in the
## beginning of this file.
#SBATCH --job-name=ep-pythia
#SBATCH --account=project_2003583
# partition explained here: https://docs.csc.fi/computing/running/batch-job-partitions/
# test = 15min, 80tasks,   2node,  382GiB max memory, 3600GiB max storage
# small= 3days, 40tasks,   1node,  382GiB max memory, 3600GiB max storage
# large= 3days, 1040tasks, 26node, 382GiB max memory, 3600GiB max storage
#SBATCH --partition=small
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type=END #uncomment to enable mail
#SBATCH --array=1-1 #defines SLURM_ARRAY_TASK_ID
#SBATCH --output=logs/output_%j.txt
#SBATCH --error=logs/errors_%j.txt

if [ "$1" == "help" ]
then
    echo "Usage: `basename $0` comment nEvts"
    exit 0
fi

if [ -z "$1" ]
then
    echo "Please give a comment to make this run unique (check '`basename $0` help' for help)"
    exit 0
fi

if [ -z "$2" ]
then
    echo "Please give input dir"
    exit 0
fi

if [ -z "$3" ]
then
    echo "Please give bmin"
    exit 0
fi

if [ -z "$4" ]
then
    echo "Please give bmax"
    exit 0
fi

outputname=$1
inputdir=$2
bmin=$3
bmax=$4

date=$(date '+%Y-%m-%d')
/users/heimarry/EventPlane/run $outputname $inputdir $bmin $bmax
sleep 1
