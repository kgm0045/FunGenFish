#! /bin/bash

source /apps/profiles/modules_asax.sh.dyn

module load sra
module load fastqc/0.10.1
module load multiqc
module load trimmomatic/0.39
module load hisat2/2.2.0
module load stringtie/2.2.1
module load gcc/9.4.0
module load python/3.10.8-zimemtc
module load samtools
module load bcftools
module load gffread

MyID=aubclsc0311          ## Example: MyID=aubclsc0311
WD=/scratch/$MyID/FinalProject             
CD=$WD/CleanData            				
REFD=$WD/Paralichthys_olivaceus_RefGenome          
MAPD=$WD/Map_STAR100           			
COUNTSD=/$WD/STAR100_Counts      
RESULTSD=/home/$MyID/FinalProject/Counts_STAR100     
REF=GCF_001970005.1_Flounder_ref_guided_V1.0_genomic.fna                  
REFA=genomic.gtf
INDEX_DIR=$WD/STAR100Index

mkdir -p $MAPD
mkdir -p $COUNTSD
mkdir -p $RESULTSD
mkdir -p $INDEX_DIR

module load star/2.7.6a
module load zlib/1.2.7 

cd ${WD}

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${INDEX_DIR} \
--genomeFastaFiles ${REFD}/${REF} \
--sjdbGTFfile ${REFD}/${REFA} \
--sjdbOverhang 100 \
--genomeSAindexNbases 13 

## Move to the directory for mapping
cd ${MAPD}

## move the list of unique ids from the original files to map
cp ${CD}/list  . 

# process the samples in the list, one by one using a while loop
while read i;
do
	STAR --genomeDir ${INDEX_DIR} \
--runThreadN 6 \
--readFilesIn  "${CD}"/"$i"_1_paired.fastq  "${CD}"/"$i"_2_paired.fastq \
--outFileNamePrefix  "${MAPD}"/"$i" \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \

mkdir "${COUNTSD}"/"$i"
	stringtie -p 6 -e -B -G  "${REFD}"/"${REFA}" -o "${COUNTSD}"/"$i"/"$i".gtf -l "$i"   "${MAPD}"/"$i"Aligned.sortedByCoord.out.bam

done<list


### The prepDE.py is a python script that converts the files in your ballgown folder to a count matrix that can be used for other programs like DESeq2. 
 ## Move to the counts directory
cd ${COUNTSD}
 ## run the python script prepDE.phy to prepare you data for downstream analysis.
cp /home/${MyID}/class_shared/prepDE.py3 .

 prepDE.py3 -i ${COUNTSD}

### copy the final results files (the count matricies that are .cvs) to your home directory. 
cp *.csv ${RESULTSD}