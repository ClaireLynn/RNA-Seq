# RNA-Seq Read Processing Guide

This detailed guide was made to instruct keen beginners to process their own RNASeq data for downstream analysis.
The bash script which generates all the directories and bash scripts for the pipline can be found in this folder [rnasetup.sh](rnasetup.sh).

Here, we will describe how to process your raw fastq files and get gene level counts for differential expression analysis.

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
mkdir claire_rna2016 

```
2) Next, find out where your files are!
They will be deposited in the ~/groupso/SO_DATA directory, then in the directory corresponding to the sequencing date.
As they are RNA-Seq samples they will begin with "R", followed by the unique 5 digit sample code, your intials eg. "CL" and end in ".fastq.gz". For example: R00197CL_R1.fastq.gz

You can use a wildcard (\*) to replace variable parts of the filename like this: R\*CL\*.fastq.gz
Check you can find your files using ls:

```
ls ~/groupso/SO_DATA/2016-01/R*CL*.fastq.gz
```

3) cd to your project and make symlinks of your fastq files.
Make a file with all symlink commands called ln.sh, then run it with bash:

```
cd claire_rna2016 

ls ~/groupso/SO_DATA/2016-01/R*CL.fastq.gz | column -t | awk '{print "ln -s "$1}' > ln.sh

bash ln.sh
```
If you make a mistake here, use ```find -type l -delete``` to delete all symlinks


4) Make prefix.txt containing the first 8 characters of the file names as follows:

```
ls *R1.fastq.gz | cut -c1-8 > prefix.txt
```

5) Run [**rnasetup.sh**](chipsetup.sh) with bash, this will make all directories and bash scripts for you to use in this guide.

You will need to give the script the genome file to map to (-g), the gtf file (-t), and the organism (-o).
Make sure the genome and gtf files exist.

For **mouse** this may be:

```
bash rnasetup.sh -g ~/groupso/REFERENCES/Mus_GRCm38.fa -t ~/groupso/REFERENCES/Mus_musculus.GRCm38.88.gtf -o mm

```
For **human** this may be:

```
bash rnasetup.sh -g ~/groupso/REFERENCES/Homo_sapiens.GRCh38.fa  -t ~/groupso/REFERENCES/Homo_sapians.GRCh38.88.gtf -o hs

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

3) cd to map, qsub **tophat2.sh** 

```
cd ../map 

qusub tophat2.sh

```


