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
SAMPLES=($InputPath/*.fastq.gz)  # Define SAMPLES as an array in Bash by using parentheses 

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}  # Get sample aligning with array index (dynamically assigned by HPC). Use ${} when the variable is part of a larger string or when there is potential ambiguity about where the variable name ends.

# set the dir to store output files
OUTDIR="/xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/fastp"  # use $HOME for passing the home directory (don't use ~, which may not be correctly expanded when passing the path to the program using it). 
mkdir -p "$OUTDIR"


# Process only Read 1 files in the SAMPLES array
if [[ $SAMPLE == *_1.fastq.gz ]]; then   # In Bash, there must be a space after [[ and before ]]
    read1=$SAMPLE
    read2=${SAMPLE/_1.fastq.gz/_2.fastq.gz}  #string substitution in Bash: ${variable/pattern/replacement}
    
    # Extract base names for read1 and read2
    base_read=$(basename "$read1" _1.fastq.gz)  # The basename command 1) takes a full file path and returns only the filename portion, and 2) removes the file suffix if provided (e.g., _1.fastq.gz). The $(command) syntax is command substitution in Bash, which captures the result of the command run inside ()

    # Run fastp program
    date
    /groups/guangyao/416a516a/modules/fastp \
          -i $read1 -I $read2 \
          -o $OUTDIR/${base_read}_1.trimmed.fastq.gz \
          -O $OUTDIR/${base_read}_2.trimmed.fastq.gz \
          -h $OUTDIR/$base_read.fastp.html \
          -j $OUTDIR/$base_read.fastp.json \
         --detect_adapter_for_pe -5 -r
          # ${var} syntax is used here to avoid the ambiguity/error of $base_read_1 being interpreted as the variable name 
    date
else
    echo "Skipped $SAMPLE: it's not an input read1 file"
fi # to end the if-then-else(optional) block in Bash
