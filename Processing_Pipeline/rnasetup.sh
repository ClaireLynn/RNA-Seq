#!/bin/bash
#$ -cwd

# get the number of samples to process
number=`wc -l prefix.txt | cut -d ' ' -f1`

# get the geneome and organism options
while getopts g:t:o: option
do
case "${option}"
in
g) genome=${OPTARG};;
t) gtf=${OPTARG};;
o) org=${OPTARG};;

esac
done

mkdir fastqc
mkdir trim
mkdir map
mkdir count

echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -N fastqc

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

fastqc $myFile\_R1.fastq.gz
fastqc $myFile\_R2.fastq.gz

exit 0;' > fastqc/fastqc.sh

echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid fastqc
#$ -N trimgalore

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

trim_galore --fastqc --paired --retain_unpaired --quality 20 \
-stringency 1  --gzip -e 0.05  \
../fastq/$myFile\_R1.fastq.gz ../fastq/$myFile\_R2.fastq.gz

exit 0;' > trim/trimgalore.sh


echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -pe smp 8
#$ -hold_jid trimgalore
#$ -N tophat2


myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

gtf='$gtf'
genome='$genome'

tophat2 -o $myFile  -p 8 --read-mismatches 2 --read-edit-dist 2 \
	--library-type fr-firststrand  -G $gtf $genome \
    ../trim/$myFile\_R1_val_1.fq.gz \
    ../trim/$myFile\_R2_val_2.fq.gz > \
    $myFile\_tophat2.log 2> $myFile\_tophat2.log.err

samtools view -f 0x02 -b $myFile/accepted_hits.bam \
| samtools sort - -n -O bam -T $myFile.nsort -o $myFile.nsort.bam


exit 0;' > map/tophat2.sh

echo \
'#!/bin/bash
#$ -cwd
#$ -t 1-'$number'
#$ -V
#$ -hold_jid tophat2
#$ -N htseq

myFile=`cat ../prefix.txt`
myFile=`head -n $SGE_TASK_ID ../prefix.txt | tail -1`
#So, the slurm ID corresponds to that line in the files

gtf='$gtf'

htseq-count --mode=union \
		-t exon \
		-i gene_id \
		-f bam \
		-r name \
		-s reverse \
		-a 30 \
		../map/$myFile.nsort.bam $gtf \
		>> $myFile.union.counts.txt


exit 0;' > count/htseq.sh
