#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name="slacker"
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=demulti_%j.out
#SBATCH --error=demulti_%j.err
#SBATCH --time=20:00:00
#SBATCH --mail-user=nszczepa@uoregon.edu
#SBATCH --mail-type=ALL 

conda activate bgmp_py37

/usr/bin/time -v ./demultiplexing.py \
-r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-i /projects/bgmp/shared/2017_sequencing/indexes.txt
