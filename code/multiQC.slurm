#!/bin/bash
 
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
### Optional. Set the job name
#SBATCH --job-name=hello_world

### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out

### REQUIRED. Specify the PI group for this job
#SBATCH --account=416a516a

### Optional. Request email when job begins and ends (if unspecified below, sent to your default university email)
### SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
### SBATCH --mail-user=<your email address>

### REQUIRED. Set the partition for your job.
#SBATCH --partition=standard

### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=1

### REQUIRED. Set the number of nodes
#SBATCH --nodes=1

### REQUIRED. Set the memory required for this job.
#SBATCH --mem=1gb

### REQUIRED. Specify the time required for this job, hh:mm:ss
#SBATCH --time=00:010:00
 
 
# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------

module load contrib

module load bjoyce3/sarawillis/multiqc/1.20

multiqc ./