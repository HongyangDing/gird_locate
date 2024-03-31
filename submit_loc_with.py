#!/bin/sh

# Lines begin with "#SBATCH" set sbatch parameters.
# Lines begin with "#" except "#!" and "#SBATCH" are comments.
# Sbatch parameters must appear before shell command. 

# Useage: sbatch intel.sh
# Output: slurm-<JOB_ID>.out

#SBATCH --get-user-env
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ddhy@pku.edu.cn

######### set job's name
#SBATCH -J with

#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

######### set TASK values(CPU CORES)
#SBATCH --ntasks 120
##SBATCH --nodes=52
##SBATCH --cpus-per-task=24
##SBATCH --exclude=cu18,cu16,cu[74-76],cu[71-72]
######### set Parallel Environment
## load environment before submitting this job
##     module load intel/2020.4.304

echo "JOB_NODELIST: ${SLURM_JOB_NODELIST}"
echo "PATH = $PATH"
echo "LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
WORKPATH=`pwd`
echo "Current Directory = $WORKPATH"

#########  execute PROGRAM_NAME
echo "Computing is started at $(date)."

$MPI_HOME/bin/mpiexec -wdir $WORKPATH python $WORKPATH/abs_loc_newnew.py
exit_code=$?

echo "Computing is stopped at $(date)."
exit $exit_code
