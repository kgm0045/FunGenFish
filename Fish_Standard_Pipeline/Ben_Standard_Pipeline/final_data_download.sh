#! /bin/bash

##############################################
## 		queue: class
##		core: 1
##		time limit (HH:MM:SS): 04:00:00 
##		Memory: 4gb
##		run on dmc
###############################################


########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsc####         ## Example: MyID=aubclsc0311

## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/$MyID/FinalProject/RawData 			
WD=/scratch/$MyID/FinalProject				
RDQ=/scratch/$MyID/FinalProject/RawDataQuality
 
##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p ${DD}
## move to the Data Directory
cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
	## -F 	Defline contains only original sequence name.
	## -I 	Append read id after spot id as 'accession.spot.readid' on defline.
	## splits the files into R1 and R2 (forward reads, reverse reads)
vdb-config --interactive

fastq-dump -F --split-files SRR13308380
fastq-dump -F --split-files SRR13308381
fastq-dump -F --split-files SRR13308382
fastq-dump -F --split-files SRR13308383
fastq-dump -F --split-files SRR13308392

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd ${WD}/${RDQ}
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*

# When finished use scpt o bring the tarballed .gz results file to your computer 
##Go to the directory on your computer that you want to copy the files to then type


