#!/bin/bash
#job name
#SBATCH --job-name=SA_resource_allocation_epi_model_calibration_runs_HIV_TB
#SBATCH --account=icrc
#SBATCH --partition=compute

## Resources
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=6GB

## You can get email notifications when your job starts and stops - useful fr long running jobs
#SBATCH --mail-user=cgreene3@uw.edu
#SBATCH --mail-type=ALL
##SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
##SBATCH --error=job.%J.err # tell it to store the error messages from the program (if it doesn't write them to normal console output) to a file called job.<assigned job muber>.err

## set array
#SBATCH --array=1-854 #change based on the number of accepted parameter (1708) sets divided by two

module load apptainer

apptainer run --bind /gscratch/icrc/cgreene3/SA_resource_allocation /gscratch/icrc/cgreene3/SA_resource_allocation/program_runs_image.sif $SLURM_ARRAY_TASK_ID

exit 0

