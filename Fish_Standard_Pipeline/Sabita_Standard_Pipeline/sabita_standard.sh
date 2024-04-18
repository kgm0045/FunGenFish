source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

WD=/scratch/aubclsc0331/Final_Project
DD=/scratch/aubclsc0331/Final_Project/RawData
RDQ=RawDataQuality
 
vdb-config --interactive

fastq-dump -F --split-files SRR13308383
fastq-dump -F --split-files SRR13308382
fastq-dump -F --split-files SRR13308381
fastq-dump -F --split-files SRR13308380
fastq-dump -F --split-files SRR13308395
fastq-dump -F --split-files SRR13308394
fastq-dump -F --split-files SRR13308393
fastq-dump -F --split-files SRR13308392


fastqc *.fastq --outdir=/scratch/aubclsc0331/Final_Project/RawDataQuality

cd /scratch/aubclsc0331/Final_Project/RawDataQuality
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*

#2. CLEAN

#! /bin/bash

######## FunGen Course Instructions ############
## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the sequence read data with Trimmomatic,
##       and then use FASTQC to evaluate the quality of the data: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
## Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: Trimmed R1 & R2 paired and unpaired reads (FASTQ)       
## FASTQC output is a folder for each file. The last line of this script will make a tarball of the output directory to bring back to your computer
##              Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
##                              Downloaded read files, R1 and R2 files for each sample if paired-end data (FASTQ)
##              Output: is a folder for each file that contains a .html file to visualize the quality, and .txt files of quality statistics.
##                      The last line of this script will make a tarball of the output directory to bring back to your computer
## For running the script on the Alabama Super Computer.
                ##For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
        ## After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
        ## then run the script by using "run_script [script name]"
        ## suggested paramenters are below to submit this script.
                ## queue: class  
                ## core: 6
                ## time limit (HH:MM:SS): 02:00:00  (may need to increase, if so run on medium queue)
                ## Memory: 12gb
                ## 
###############################################

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

# Modules
        #  load the module
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

## STOP. You need to replace the [number] with YOUR paths to 
##       make variables for your ASC ID so the directories are automatically made in YOUR directory
MyID=aubclsc0331                        ## Example: MyID=aubtss

# Variables: raw data directory (DD), working directory(WD), Quality after cleaning (PCQ), name of file containing the adpaters.
WD=/scratch/$MyID/FinalProject                         ## Example: WD=/scratch/$MyID/PracticeRNAseq
DD=/scratch/$MyID/FinalProject/RawData                          ## Example: DD=/scratch/$MyID/PracticeRNAseq/RawData
CD=/scratch/$MyID/FinalProject/CleanData                          ## Example: CD=/scratch/$MyID/PracticeRNAseq/CleanData
PCQ=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing.

                               ## In the future, for your data, you will likely need to edit this for other projects based on how your libraries 
                                ## were made to search for the correct adapters for your project
                                
## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir ${CD}
mkdir ${WD}/${PCQ}

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
        ## CHECK: You may need to edit this path for the file that is in the class_shared directory from your account.
cp /home/${MyID}/class_shared/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code


while read i
do

       ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
	       java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
					PE -threads 6 -phred33 \
        	"$i"_1.fastq "$i"_2.fastq  \
       	 ${CD}/"$i"_1_paired.fastq ${CD}/"$i"_1_unpaired.fastq  ${CD}/"$i"_2_paired.fastq ${CD}/"$i"_2_unpaired.fastq \
       	 ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.


        ############## FASTQC to assess quality of the Cleaned sequence data
        ## FastQC: run on each of the data files that have 'All' to check the quality of the data
        ## The output from this analysis is a folder of results and a zipped file of results

fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}

done<list                       # This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd ${WD}/${PCQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*

## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate the quality of the data.


#3. MAPING

#!/bin/sh
 
######### FunGen Course Instructions ############
## Purpose: The purpose of this script is to 
##    Use HiSat2 to index your reference genome and then map your cleaned (paired) reads to the indexed reference
##              First need to use gffread to convert annotation file from .gff3 to .gft formate
##              Use Stringtie to count the reads mapped to genes and transcripts, defined in this case by the genome annotation file
##              use the python script to take the Stringtie results to make two counts matricies, one at the gene level and one at the transcript l$
## HiSat2  Indexing   InPut: Reference genome file (.fasta), and annotation file (.gff3) (Optional)
##                    Output: Indexed genome 
## HiSat2 Mapping     Input: Cleaned read files, paired (.fasq); Indexed genome
##                    Output: Alignment .sam files  
## Samtools  Convert .sam to .bam and sort          Input: Alignment files,  .sam
##                                                  Output: Sorted  .bam files
## Stringtie  Counting reads  Input: sorted .bam file
##                            Output:  Directories of counts files for Ballgown (R program for DGE)
##              prepDE.py    Python script to create a counts matrics from the Stringtie output.  Inputs: Directory from Stringtie
##                                                                                                Output:  .csv files of counts matrix
## For running the script on the Alabama Super Computer.
##  For more information: https://hpcdocs.asc.edu/content/slurm-queue-system
##  After you have this script in your home directory and you have made it executable using  "chmod +x [script name]", 
##  then run the script by using "run_script [script name]"
##  suggested paramenters are below to submit this script.
##    queue: class or medium
##    core: 6
##    time limit (HH:MM:SS): 04:00:00 
##    Memory: 12gb
##    
###############################################

#### Load all the programs you are going to use in this script.
source /apps/profiles/modules_asax.sh.dyn
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread
#module load gffcompare


#  Set the stack size to unlimited
ulimit -s unlimited

# Turn echo on so all commands are echoed in the output log
set -x

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
  ## Replace the [#] with paths to define these variable
  MyID=aubclsc0331

WD=/scratch/$MyID/FinalProject           ## Example:/scratch/$MyID/PracticeRNAseq  
CD=/scratch/$MyID/FinalProject/CleanData            ## Example:/scratch/$MyID/PracticeRNAseq/CleanData   #   *** This is where the cleaned paired f$
REFD=/scratch/$MyID/FinalProject/ReferenceGenome          ## Example:/scratch/$MyID/PracticeRNAseq/DaphniaRefGenome    # this directory contains th$
MAPD=/scratch/$MyID/FinalProject/Map_HiSat2           ## Example:/scratch/$MyID/PracticeRNAseq/Map_HiSat2      #
COUNTSD=/scratch/$MyID/FinalProject/Counts_StringTie	   ## Example:/scratch/$MyID/PracticeRNAseq/Counts_StringTie
RESULTSD=/home/$MyID/FinalProject/Counts_H_S_2024      ## Example:/home/aubtss/PracticeRNAseq/Counts_H_S

               ## This is what the "easy name" will be for the genome reference

## Make the directories and all subdirectories defined by the variables above
mkdir -p $REFD
mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD

##################  Prepare the Reference Index for mapping with HiSat2   #############################
cd $REFD


###  Identify exons and splice sites on the reference genome
gffread genome.gff -T -o Paralichthys_olivaceus.gtf              ## gffread converts the annotation file from .gff3 to .gft formate$
hisat2_extract_splice_sites.py Paralichthys_olivaceus.gtf > Paralichthys_olivaceus.ss
hisat2_extract_exons.py Paralichthys_olivaceus.gtf> Paralichthys_olivaceus.exon

#### Create a HISAT2 index for the reference genome. NOTE every mapping program will need to build a its own index.
hisat2-build --ss Paralichthys_olivaceus.ss --exon Paralichthys_olivaceus.exon Paralichthys_olivaceus.fasta Paralichthys_index

########################  Map and Count the Data using HiSAT2 and StringTie  ########################

# Move to the data directory
cd ${CD}  #### This is where our clean paired reads are located.

## Create list of fastq files to map.    Example file format of your cleaned reads file names: SRR629651_1_paired.fastq SRR629651_2_paired.fastq
## grab all fastq files, cut on the underscore, use only the first of the cuts, sort, use unique put in list
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list Example: SRR629651

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
mv ${CD}/list  . 

## process the samples in the list, one by one using a while loop
while read i;
do
  ## HiSat2 is the mapping program
  ##  -p indicates number of processors, --dta reports alignments for StringTie --rf is the read orientation
   hisat2 -p 6 --dta --phred33       \
    -x "${REFD}"/Paralichthys_index	 \
    -1 "${CD}"/"$i"_1_paired.fastq  -2 "${CD}"/"$i"_2_paired.fastq	\
    -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: SRR629651.sam; Output: SRR629651_sorted.bam
  samtools view -@ 6 -bS "$i".sam > "$i".bam  

    ###  This is sorting the bam, using 6 threads, and producing a .bam file that includes the word 'sorted' in the name
  samtools sort -@ 6  "$i".bam  -o  "$i"_sorted.bam

    ### Index the BAM and get mapping statistics, and put them in a text file for us to look at.
  samtools flagstat   "$i"_sorted.bam   > "$i"_Stats.txt

### Stringtie is the program that counts the reads that are mapped to each gene, exon, transcript model. 
  ### The output from StringTie are counts folders in a directory that is ready to bring into the R program Ballgown to 
  ### Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  ### eAB options: This will run stringtie once and  ONLY use the Ref annotation for counting readsto genes and exons 
  
mkdir "${COUNTSD}"/"$i"
stringtie -p 6 -e -B -G  "${REFD}"/Paralichthys_olivaceus.gtf -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"_sorted.bam

done<list

#####################  Copy Results to home Directory.  These will be the files you want to bring back to your computer.
### these are your stats files from Samtools
cp *.txt ${RESULTSD}

### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix. 
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv ${RESULTSD}

## move these results files to your personal computer for downstream statistical analyses in R studio.
