#!/bin/bash

### Don't change unless you know why
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
### Each task requires one cpu by default; include cpus-per-task if more are required
#SBATCH --mem=10gb
#SBATCH --time=02:00:00

### Change as you see fit
#SBATCH --account=416a516a
#SBATCH --mail-type=END

### Featurecount job part
courseData="/groups/guangyao/416a516a"
MP="$courseData/modules/subread-2.0.3-Linux-x86_64/bin"
annotation="/groups/guangyao/416a516a/mouse_genome.annotation/gencode.vM36.primary_assembly.annotation.gtf"
IP=$(ls /xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/STAR/*.bam) # The input files for featureCounts should be the output files (mapped reads) of STAR mapping, the previous step in the pipeline

### set output file dir
OP="/xdisk/guangyao/416a516a/516_projects/Group_5/bryanjacobs/featureCounts" ## replace OP with your own output dir
mkdir -p $OP

date
$MP/featureCounts \
-T 1 \
-p --countReadPairs \
-a $annotation \
-o $OP/featureCounts_PE.txt \
$IP
date