#!/bin/bash

### SLURM directives. Must be place before any shell commands or variable definitions to be recognized
# Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=5gb
#SBATCH --time=00:20:00

# Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END

# Define the job array
#SBATCH --array=0-105  # Adjust range based on number of samples: <0 - (number of files -1)>

### Bash commands
# Define input files and directories
InputPath="/xdisk/guangyao/416a516a/516_projects/Group_5/run_files"
SAMPLES=($InputPath/*.fastq.gz)  # # Define SAMPLES as an array in Bash by using parentheses 

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}  # Get sample aligning with job array index (dynamically assigned by HPC). Use ${} when the variable is part of a larger string or when there is potential ambiguity about where the variable name ends.

# set the dir to store output files
OUTDIR="/xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/fastQC" # use $HOME for passing the home directory (don't use ~, which may not be correctly expanded when passing the path to the program using it).  

module load fastqc	### Load the HPC-installed module/program

fastqc $SAMPLE -o $OUTDIR

# -o $OUTDIR is important: fastqc otherwise generate the output files in the input file dir by default