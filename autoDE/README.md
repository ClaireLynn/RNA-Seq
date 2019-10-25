
# autoDE

R scripts which accept arguments and perform differential expression analysis without the need for writing any code.
These scripts are perfect for when there is a simple model, or want to try many different comparisons and you need fast results!


Run like so: 
```
Rscript --vanilla autoDE.R -n my_analysis -s samples.txt -m ~condition+type -c contrasts.txt -o mouse
```
 
To get information on the input required use the --help / -h option:
```
Rscript --vanilla autoDE.R --help 
```

You must install the following packages from bioconductor (bioconductor.org, BiocManager::install()): 
    DESeq2
    TxDb.Mmusculus.UCSC.mm10.knownGene
    biomaRt 
    org.Mm.eg.db
From CRAN (install.packages())
RColorBrewer
ggplot2 
ggtern
gplots
grid

This script will spit out the following results files:
name_PCA.pdf-This file is a PCA plot of your samples 
contrast,.csv -This is your results, with normalised counts, ensembl and gene symbols
contrast,_volc_plot.pdf -Volcano plot for particular contrast, threshold passing genes labelled
 
If organism is mouse, you will also find the following additional files:
 
contrast_human_symbols.csv -Results With human symbols 
contrast_human_symbols_nodup.csv -Results with human symbols, duplicates human genes removed 
contrast_human_symbols.rnk -human symbols, hgnc symbol and logfc ready for GSEA!
