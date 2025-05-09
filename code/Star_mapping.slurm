#!/bin/bash

### Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --mem=32gb

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
### Each task requires one cpu by default; include cpus-per-task if more are required

#SBATCH --time=48:00:00

### Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END

### STAR mapping job part

# Define annotation file (used for downstream analysis like featureCounts)
annotation="/groups/guangyao/416a516a/mouse_genome.annotation/gencode.vM36.primary_assembly.annotation.gtf"

### set input file dir
courseData="/groups/guangyao/416a516a"

InputDir="/xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/fastp"
OUTDIR="/xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/STAR"
mkdir -p "$OUTDIR"

# Define the list of input files (Read 1 files, automatically handled by STAR for Read 2)
SAMPLES=($InputDir/*_1.trimmed.fastq.gz)

module load star

# Loop through each Read 1 file and map with Read 2
for read1 in "${SAMPLES[@]}"; do
    # Get the base name of the Read 1 file
    base_read=$(basename "$read1" _1.trimmed.fastq.gz)
    
    # Define the corresponding Read 2 file
    read2="${InputDir}/${base_read}_2.trimmed.fastq.gz"
    
    # Run STAR for each pair of Read 1 and Read 2
    date
    STAR \
    --runThreadN 2 \
    --genomeDir /groups/guangyao/416a516a/mouseIndex \
    --readFilesIn $read1 $read2 \
    --readFilesCommand zcat \
    --outFileNamePrefix $OUTDIR/${base_read}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif
    date
done