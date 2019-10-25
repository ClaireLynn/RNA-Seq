# RNA-Seq Read Processing Guide

This detailed guide was made to instruct keen beginners to process their own RNASeq data for downstream analysis.
The bash script which generates all the directories and bash scripts for the pipline can be found in this folder: **rnasetup.sh**.

Here, we will describe how to process your raw fastq files and get peaks and track files for IGV.

## Notes on bash scripts

First, a note on bash scripts using the Sun Grid Engine (SGE) environment. 
At the top of every script are the qsub options:

```
#!/bin/bash
#$ -cwd
#$ -t 1-12
#$ -V
#$ -hold_jid fastqc
#$ -N trimgalore
#$ -l h_rt=10:00:00

```

```#!/bin/bash``` allows this to be submitted as a batch job using ```qsub```.
```-cwd``` enforces that the script is run in the current working directory.
```-t``` is the task IDs to run, we use a prefix file to cycle through these IDs creating a job array in this part of each script :

```
myFile=`cat prefix.txt`
myFile=`head -n $SGE_TASK_ID prefix.txt | tail -1`
```

```-V``` allows bash to find programs in your path
```-hold_jid``` allows you to hold the job in the queue until the specified job has completed
```-N``` allows you to name the job

If your jobs are queued but not running for a long time try adding ```-l h_rt=06:00:00``` after qsub for a max run time of 6 hours,
the default option reserves a run time of 10 hours.

You may also add additional options when you run ```qsub```, reference: http://gridscheduler.sourceforge.net/htmlman/manuals.html for more options.



## Preparation

1) In your scratch directory, make a directory for your project

``` 
mkdir claire_chip2016 

```
2) Next, find out where your files are!
They will be deposited in the ~/groupso/SO_DATA directory, then in the directory corresponding to the sequencing date.
As they are ChIP samples they will begin with "C", followed by the unique 5 digit sample code, your intials eg. "CL" and end in ".fastq.gz". For example: C00447CL_R1.fastq.gz

You can use a wildcard (\*) to replace variable parts of the filename like this: C\*CL\*.fastq.gz
Check you can find your files using ls:

```
ls ~/groupso/SO_DATA/2016-04/C*CL*.fastq.gz
```

3) cd to your project and make symlinks of your fastq files.
Make a file with all symlink commands called ln.sh, then run it with bash:

```
cd claire_chip2016 

ls ~/groupso/SO_DATA/2018-03/C*CL.fastq.gz | column -t | awk '{print "ln -s "$1}' > ln.sh

bash ln.sh
```
If you make a mistake here, use ```find -type l -delete``` to delete all symlinks


4) Make prefix.txt containing the first 8 characters of the file names as follows:

```
ls *R1.fastq.gz | cut -c1-8 > prefix.txt
```

5) Run **chipsetup.sh** with bash, this will make all directories and bash scripts for you to use in this guide.

You will need to give the script the genome file to map to (-g) and the organism (-o).
Make sure the genome file exists.

For **mouse** this may be:

```
bash chipsetup.sh -g ~/groupso/REFERENCES/Mus_GRCm38.fa -o mm

```
For **human** this may be:

```
bash chipsetup.sh -g ~/groupso/REFERENCES/Homo_sapiens.GRCh38.fa -o hs

```

## Processing

Now we begin actually processing these data. Each job must be ran in the correct order but may be submitted before previous jobs have ended.

1) cd to fastqc and run **fastqc.sh** with qsub.

```
cd fastqc 

qsub fastqc.sh

```

2) cd to trim, qsub **trimgalore.sh**

```
cd ../trim 

qusub trimgalore.sh

```

3) cd to map, qsub **bowtie2.sh** 

```
cd ../map 

qusub bowtie2.sh

```

4) Remaining in the map directory, qsub **rmdup.sh** to remove duplicates

``` 
qusub rmdup.sh

```

## Peak calling and visualisation

1) cd to peakcall and make two files:

**"prefix_test.txt" - contains all prefixes for ChIP samples.**

**"prefix_ctrl.txt" - contains all prefixes for input controls (or igg)**, these must be in the same order corresponding to the chip samples. 
eg: If C00005CL is the input control for ChIP samples: C00001CL/C00002CL and C00006CL is the input control for C00003CL/C00004CL the contents of the files would be as follows:

        prefix_test.txt:                prefix_ctrl.txt:
        C00001CL                        C00005CL
        C00002CL                        C00005CL
        C00003CL                        C00006CL
        C00004CL                        C00006CL

2) qsub **macs2.sh** for narrow peak calling (best for TFs) or **macs2_broad.sh** for broad peak calling (best for histone marks).

```
qsub masc2.sh
```

3) To make bigwig files for viewing in IGV, you require a file inside the peakcall directory called **"chrom.sizes.txt"** which contains the chromosome sizes of your organism. You can find this online (make sure the chromosome names match the reference files!) or create it from your bam
file header! Extract with samtools view and then use sed to remove the
unwanted characters e.g.:

```
samtools view -H ../map/C00001CL.bam | head -n -1 | sed 's/@SQ\t//g; s/SN://g; s/\t/ /g; s/LN://g' | tail -n +2 > chrom.sizes.txt
```
4) Then you can run **makebigwigs.sh** remaining in the peakcall directory

```
qsub makebigwigs.sh

```

5) *Optional*: Use bedtools to look for differences between peak files. You need to
decide which are the best options to use as this is project specific.

http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
eg.
```
bedtools intersect -wa -a C00012RK_peaks.broadPeak -b C00016RK_peaks.broadPeak > wa.12v16.broadPeak
```

http://bedtools.readthedocs.io/en/latest/content/tools/subtract.html
eg.
```
bedtools subtract -A -a C00012RK_peaks.broadPeak -b C00016RK_peaks.broadPeak > 12v16.subtract.broadPeak
```


6) Any further steps will now be completed in R using bioconductor packages. This is highly project specific.
       
