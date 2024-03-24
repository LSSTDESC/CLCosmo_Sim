#!/bin/bash

# SLURM options:

#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --output=serial_test_%j.log   # Standard output and error log

#SBATCH --partition=htc               # Partition choice (htc by default)

#SBATCH --ntasks=1                    # Run a single task
#SBATCH --mem=2000                    # Memory in MB per default
#SBATCH --time=1-00:00:00             # Max time limit = 7 days

#SBATCH --mail-user=<email address>   # Where to send mail
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

#SBATCH --licenses=sps                # Declaration of storage and/or software resources

# Commands to be submitted: