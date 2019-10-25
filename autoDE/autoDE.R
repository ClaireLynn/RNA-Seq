#!mnt/lustre/groups/groupso/bin/ Rscript
cat("
      
                                                                                                                     
                   autoDE was written Claire Lynn, happy to help: clairepaulalynn@gmail.com                            
                                                                                                                     
Run like so:                                                                                               
   Rscript --vanilla autoDE.R -n my_analysis -s samples.txt -m ~condition+type -c contrasts.txt -o mouse  
                                                                                                                     
   To get information on the input required use the --help / -h option:                                      
      Rscript --vanilla autoDE.R --help                                                                     
                                                                                                                     
      You must install the following packages from bioconductor 
      (bioconductor.org, BiocManager::install()):     
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
        name_PCA.pdf               -This file is a PCA plot of your samples                                 
        contrast,.csv              -This is your results, with normalised counts, ensembl and gene symbols  
        contrast,_volc_plot.pdf    -Volcano plot for particular contrast, threshold passing genes labelled  
                                                                                                                     
      If organism is mouse, you will also find the following additional files:                                    
                                                                                                                     
        contrast_human_symbols.csv       -Results With human symbols                                               
        contrast_human_symbols_nodup.csv -Results with human symbols, duplicates human genes removed               
        contrast_human_symbols.rnk       -human symbols, hgnc symbol and logfc ready for GSEA!                     
                                                                                                                     
                                                                                                                     
")
      
suppressMessages(library("optparse"))

# Specify options for this script: name, samples, model, contrasts

option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="The name of your analysis, output will be named by this",
              metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL, 
              help="Samples file name, contains tab separated columns with data file names and metadata for your model, e.g.:
                                 
                                type  replicate                  
                  countfile1.txt big   1
                  countfile2.txt small 1
                  countfile3.txt big   2
                  countfile4.txt small 2

              >>> The row names MUST be the file names <<<",
		          metavar="character"),
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="The DESeq2 model to be used for the differential analysis as a formula: ~condition1+condition2 for more information about 
              the formula please refer to the DESeq2 documentation http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html",
		          metavar="character"),
  make_option(c("-c", "--contrasts"), type="character", default=NULL, 
              help="give a text file with a line separated list of contrasts in the following format: condition,groupA,groupB
              The condition MUST match a column in --samples, groupA / groupB MUST match a group in the model",
		          metavar="character"),
  make_option(c("-o", "--organism"), type="character", default=NULL, 
              help="human or mouse",
              metavar="character")
) 

opt_parser = OptionParser(option_list=option_list)
# Give the users arguments to the object "opt"
opt = parse_args(opt_parser)



# Print error messages if arguments are missing

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Must supply a name!", call.=FALSE)
}

if (is.null(opt$samples)){
  print_help(opt_parser)
  stop("Must supply samples file!", call.=FALSE)
}

if (is.null(opt$model)){
  print_help(opt_parser)
  stop("Must supply a model!", call.=FALSE)
}

if (is.null(opt$contrasts)){
  print_help(opt_parser)
  stop("Must supply contrasts!", call.=FALSE)
}

if (is.null(opt$organism)){
  print_help(opt_parser)
  stop("Must supply organism!", call.=FALSE)
}

# Get the information from the user's inputs 
samples <- read.table(opt$samples, header=T, stringsAsFactors = TRUE)
print(samples)

#convert any numbers in this table to factors
samples[sapply(samples, is.integer)] <- lapply(samples[sapply(samples, is.integer)], as.factor)

files <- row.names(samples)

model<-opt$model
name <-opt$name
org <-opt$organism

contrasts<-scan(opt$contrasts, what="",sep="\n")

# Read in the count files, assign them to objects named by their file name
for (file in files)
{
assign(file,read.table(file)[,2])
}

# Get the list of objects generated as a result of the above code (read count file objects)
objectlist<- lapply(files, get)
names(objectlist)<- files

# Use the object list to column bind the counts
read.count.table<- do.call(cbind,objectlist)

# Get the gene names
file<-files[1]
gene.names = as.character((read.table(file)[,1]))

# Name the samples from the first 8 characters of the supplied file names (HOPEFULLY, this is the 8 chararacter sample code!!)
sample.names <- substr(files,1 , 8)
row.names(samples) <- sample.names

# Name the rows and columns by gene names and sample names respectively
dimnames(read.count.table)<-list(gene.names,sample.names)

# Filter your count files for ambiguous counts and ammend gene.names ("_" found in htseq, "+" found in MZ count files) 
read.count.table <- read.count.table[-grep("\\_|\\+",row.names(read.count.table)),]
gene.names = as.character(row.names(read.count.table))

# Load all the packages and suppressMessages
suppressMessages(library(DESeq2))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggtern"))
suppressMessages(library("ggrepel"))
suppressMessages(library("gplots"))
suppressMessages(library(grid))
suppressMessages(library("biomaRt"))
suppressMessages(library("org.Mm.eg.db"))

# Write sessionInfo - very useful for replication!
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


# Make the DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=as.matrix(read.count.table),
	 colData=samples,
	 design=as.formula(model))

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# Get the rlog transformed values (useful for PCA!)
rld <- rlogTransformation(dds)
rld2 <- rlogTransformation(dds, blind=FALSE)
write.table(assay(rld),paste0(name,"_rld.txt"))
write.table(assay(rld2),paste0(name,"_rld_model.txt"))

# Get the vst transformed values (useful for heatmaps)
vst <- varianceStabilizingTransformation(dds)
vst2 <- varianceStabilizingTransformation(dds, blind=FALSE)
write.table(assay(vst),paste0(name,"_vst.txt"))
write.table(assay(vst2),paste0(name,"_vst_model.txt"))

# Make PCA plot - plot_myPCA is a modified version of plotPCA so that the shape and colours of points are informed by 
# columns of the sample information and points are labelled by their sample names 

plot_myPCA <- function(rld,name,Height,Width,color,shape=NULL,samples=samples){
data <- plotPCA(rld,intgroup=c(names(samples)),returnData=TRUE,ntop=5000)
percentVar <- round(100*attr(data,"percentVar"))

plot <- ggplot(data,aes(PC1, PC2, label=substr(row.names(samples),1,8),
		color=eval(parse(text=color)), shape=eval(parse(text=shape))))+
		geom_point(size=5)+
		xlab(paste0("PC1: ", percentVar[1],"% variance"))+
		ylab(paste0("PC2: ",percentVar[2],"% variance"))+
		scale_colour_discrete(name=color)+
		scale_shape_discrete(name=shape)+
		geom_text_repel(color="black", size=3, point.padding=NA)
plot +theme_classic()
ggsave(filename=paste0(name,"_PCA.pdf"), plot = last_plot(), scale = 1, 
width = Width, height = Height, units ="cm")
}


# Get the colour and shape of the PCA plot points, if there is more than one sample information column,
# shape is informed by the second column.
color<-names(samples)[1]


if (ncol(samples) > 1) {
  shape <- names(samples)[2]
} else {
   shape <- names(samples)[1]
}

# Run the PCA function
print("Saving PCA plot as PDF")
plot_myPCA(rld=rld,name=name,Height=10,Width=12,color=color,shape=shape, samples=samples)


######### Get gene symbols and human orthologs #########




if(opt$organism=="mouse"){
# Name the biomart database you want to use
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset= "mmusculus_gene_ensembl")
# Fetch the mgi_symbols and produce a table with ensembl_gene_id | mgi_symbol
mouse_symbols <- getBM(attributes = c("ensembl_gene_id","mgi_symbol","entrezgene_id"),
                       filters='ensembl_gene_id',
                       values=gene.names,
                       mart=ensembl)
# get the human/mouse genes
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human_symbols = getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), 
                       mart = mouse, attributesL = c("hgnc_symbol", "ensembl_gene_id"), 
                       martL = human, uniqueRows=T)
# rename the columns
human_symbols <- rename(human_symbols, 
			c("Gene.stable.ID"="Ensembl.Gene.ID.Mouse",
			"Gene.stable.ID.1"="Ensembl.Gene.ID.Human"))
}

if(opt$organism=="human"){
# Name the biomart database you want to use
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset= "hsapiens_gene_ensembl")
human_symbols <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
                       filters='ensembl_gene_id',
                       values=gene.names,
                       mart=ensembl)
}

##### Perform DESeq function #####
print("Now Performing DESeq!")
dds <- DESeq(dds)

if(opt$organism=="mouse"){
norm_counts<-merge(mouse_symbols,as.data.frame(counts(dds,normalized=TRUE)),by.x=1,by.y=0, all.y=TRUE)
}

if(opt$organism=="human"){
norm_counts<-merge(human_symbols,as.data.frame(counts(dds,normalized=TRUE)),by.x=1,by.y=0, all.y=TRUE)
}


write.csv(norm_counts, file=paste0(name,"_norm_counts.csv"),row.names=FALSE,quote=FALSE)

# Get your DESeq results
# The following function extracts the results from your deseq object for a given contrast
# It will write:
# a deseq results table with the MGI symbols and normalised read counts appended (csv)
# a deseq results with the HGNC symbols and normalised read counts appended (csv)
# the above, fixed for many to many and one to many ortholog mapping (csv)
# the above, only the HGNC symbol and logfc columns, ready for GSEA (csv)
# a volcano plot labelled and unlabelled (pdf)

SaveDeseqResults <- function(name,dds,contrast,
				gene.names=gene.names,
			 	mouse=mouse_symbols,
			 	human=human_symbols){
	
	res<-results(dds, contrast=unlist(strsplit(contrast, split=",")))
	res<-merge(as.data.frame(res),as.data.frame(counts(dds, normalized=TRUE)),by.x=0,by.y=0)
	res<-rename(res, c("Row.names"="Ensembl.Gene.ID"))
	write.csv(as.data.frame(res),file=paste0(name,".csv"), row.names=FALSE,quote=FALSE)	

	res_m<-merge(as.data.frame(mouse_symbols),as.data.frame(res),by.x=1,by.y=1, all.y=TRUE)
	res_m<-na.omit(res_m)
	res_m<-res_m[order(res_m$padj),]
	write.csv(as.data.frame(res_m),file=paste0(name,"_symbols.csv"), row.names=FALSE,quote=FALSE)

	res_h<-merge(as.data.frame(human),as.data.frame(res),by.x=1,by.y=1, all.y=TRUE)
	res_h[res_h == ""] <- NA
	res_h<-na.omit(res_h)
  	res_h<-res_h[order(res_h$padj),]
	write.csv(as.data.frame(res_h),file=paste0(name,"_human_symbols.csv"), row.names=FALSE,quote=FALSE)

		########## FIXING DUPLICATES FOR GSEA #########
  # Logic:
	# many human to 1 mouse will be given the same logfc
	# many mouse to 1 human have an averaged logfc 
	# many mouse to many human have averaged logfc of all mouse genes which map to all human mappings

	# first antilog the LogFC values, then average them :
	res_h$FoldChange <- 2^res_h$log2FoldChange
	# Aggregate the foldchange values (mean) by the HGNC symbols and 
	# then merge the aggregated values with the rest of the table
	res_h <- merge(aggregate(FoldChange ~ HGNC.symbol, data = res_h, FUN = mean), res_h, all.x=TRUE)
	
	# log2 all the fold changes to get the log2fold changes: this will add the log2foldchanges for the aggregated values 
	res_h$log2FoldChange <- log2(res_h$FoldChange)

	# Select the HGNC column, followed by the log2foldchanges for a ranked list
	write.csv(res_h,file=paste0(name,"_human_symbols_dedup.csv"),row.names=FALSE,quote=FALSE)
	res_h<-res_h[,c(1,7)]
	write.table(res_h,file=paste0(name,"_human_symbols_dedup.rnk"),row.names=FALSE,
	col.names=FALSE, quote=FALSE, sep="\t", )
	

  # Make volcano plot
	res_m$threshold = as.factor(abs(res_m$log2FoldChange) > 1 & res_m$padj < 0.05)
	print(head(res_m))
	

	volc <- ggplot(data=as.data.frame(res_m), 
	aes(x=log2FoldChange, y=-log10(padj))) +
        geom_point(aes(col=threshold),alpha=0.4,size=4) +
        xlab("Log Fold Change") + ylab("-log10 Adjusted P.Value")+
        geom_text_repel(inherit.aes = FALSE, 
	                     data=head(res_m[abs(res_m$log2FoldChange)>1,],50),
        aes(x=log2FoldChange, y=-log10(padj),label = mgi_symbol),
	         size=3, segment.color ="grey", segment.alpha=0.5)
	
	volc + theme_classic(base_size = 12)
	ggsave(filename=paste0(name,"_volc_plot.pdf"),
	 plot = last_plot(), scale = 1, width = 16, height = 18, units ="cm")
	
}

SaveDeseqResultsh <- function(name,dds,contrast,
        gene.names=gene.names,
        human=human_symbols){

  res<-results(dds, contrast=unlist(strsplit(contrast, split=",")))
  res<-merge(as.data.frame(res),as.data.frame(counts(dds, normalized=TRUE)),by.x=0,by.y=0)
  res<-rename(res, c("Row.names"="Ensembl.Gene.ID"))
  write.csv(as.data.frame(res),file=paste0(name,".csv"), row.names=FALSE,quote=FALSE) 

  res_h<-merge(as.data.frame(human),as.data.frame(res),by.x=1,by.y=1, all.y=TRUE)
  res_h[res_h == ""] <- NA
  res_h<-na.omit(res_h)
  res_h<-res_h[order(res_h$padj),]
  write.csv(as.data.frame(res_h),file=paste0(name,"_human_symbols.csv"), row.names=FALSE,quote=FALSE)

  # Select the HGNC column, followed by the log2foldchanges for a ranked list
  write.csv(res_h,file=paste0(name,"_human_symbols_dedup.csv"),row.names=FALSE,quote=FALSE)
  res_h<-res_h[,c(1,7)]
  write.table(res_h,file=paste0(name,"_human_symbols_dedup.rnk"),row.names=FALSE,
  col.names=FALSE, quote=FALSE, sep="\t", )
  
  # Make volcano plot
  res_h$threshold = as.factor(abs(res_h$log2FoldChange) > 1 & res_h$padj < 0.05)
  print(head(res_h))
  
  volc <- ggplot(data=as.data.frame(res_h), 
  aes(x=log2FoldChange, y=-log10(padj))) +
        geom_point(aes(col=threshold),alpha=0.4,size=4) +
        xlab("Log Fold Change") + ylab("-log10 Adjusted P.Value")+
        geom_text_repel(inherit.aes = FALSE, 
                       data=head(res_h[abs(res_h$log2FoldChange)>1,],50),
        aes(x=log2FoldChange, y=-log10(padj),label = hgnc_symbol),
           size=3, segment.color ="grey", segment.alpha=0.5)
  
  volc + theme_classic(base_size = 12)
  ggsave(filename=paste0(name,"_volc_plot.pdf"),
   plot = last_plot(), scale = 1, width = 16, height = 18, units ="cm")
  
}

#Apply this function to how ever many contrasts you have:
if(opt$organism=="mouse"){
for (contrast in contrasts)
  {
    SaveDeseqResults(dds=dds, 
                      name=paste0(name,"_",gsub(",","_",contrast)), 
                      contrast=contrast, 
                      mouse=mouse_symbols, 
                      human=human_symbols, 
                      gene.names=gene.names)
  }
}

if(opt$organism=="human"){
for (contrast in contrasts)
  {
    SaveDeseqResultsh(dds=dds, 
                      name=paste0(name,"_",gsub(",","_",contrast)), 
                      contrast=contrast, human=human_symbols, 
                      gene.names=gene.names)
  }
}

