#---------Bulk RNAseq Analysis

#Set up the working directory----
setwd("C:\\Users\\Jehad Yasin\\OneDrive\\Desktop\\Cancer Genomics")

#Loading required packages----
#install.packages("R.utils")
library(R.utils)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rsubread")
library(Rsubread)

#install.packages("data.table")
#Mapping, quantification and variant analysis
library(data.table)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RUVSeq")
#Remove unwanted variation from RNA-Seq
library("RUVSeq")

#(if the above installation of RUVseq didn't work, try this)
#source("http://bioconductor.org/biocLite.R")
#biocLite("RUVSeq")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
#Differential GeneExp based on (-) binomial dist.
library(DESeq2)

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
#volcano plots
library(EnhancedVolcano)

#install.packages("pheatmap")
library(pheatmap)

#install.packages("RColorBrewer")
library(RColorBrewer)

#install.packages("ggplot2")
library(ggplot2)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Rqc")
#quality control tool
library(Rqc)

#Downloading FastQ files----

#-------------------------- Downloaded Manually

#Downloading Genome file----
#Reference genome (source: ensembl genome browser)
#url<-"ftp://www.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
#destination<-"Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
#download.file(url,destination)
#gunzip(destination)

#Downloading GTF file----
#GTF File containing indexes (source:ensembl genome browser)
#url<-"ftp://www.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
#destination<-"Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"
#download.file(url,destination)
#gunzip(destination)



#rqc(path = "C://Users//Jehad Yasin//OneDrive//Desktop//Cancer Genomics", pattern = ".fastq.gz")

#Building Index----
#saving index to a file named "Sc_full_index_rsubread"
buildindex("Sc_full_index_rsubread",
           "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa",
           indexSplit=F)

#check that we have all files
reads1 <- list.files(pattern = "_1.fastq.gz$" )
reads2 <- list.files(pattern = "_2.fastq.gz$" )
all.equal(length(reads1),length(reads2))

#Performing Alignment----
align(index="Sc_full_index_rsubread",
      readfile1=reads1,
      readfile2=reads2,
      input_format="gzFASTQ",
      output_format="BAM",
      nthreads=10)

#Checking the BAM files generated----
#2 BAM files produced, one for each sample. 
#Samples: SRR5924196 and SRR5924198 (e.g. Normal vs Tumor)
bam.files <- list.files(pattern = ".BAM$", full.names = TRUE)
bam.files

#Checking the mapping quality----
props <- propmapped(files=bam.files)
props

#Generating Feature counts----
#BAM files to Count files by aid of exon-GTF indexing
#fclim is a list of lists
#double indexing [[]] for list in list
fcLim <- featureCounts(files = bam.files,
                       annot.ext="Saccharomyces_cerevisiae.R64-1-1.96.gtf",
                       GTF.featureType="exon",
                       GTF.attrType="gene_id",
                       isGTFAnnotationFile=TRUE,
                       isPairedEnd = TRUE)

#try View "fc"
#fc is "counts" extract from fclim
fc <- data.frame(fcLim[["counts"]])
colnames(fc) <- c("Normal", "Tumor")

#Working on actual data----
f <- read.csv("GSE143630_RCC_htseq_counts.txt",
              sep=" ", 
              row.names = 1)

#T1 stage count---- (Stage 1)
count_of_T1=length(grep(x = colnames(f),
                            pattern = "^T1."))
count_of_T1 #--number of samples with Stage 1

#T2 stage count---- (Stage 2)
#remember "grep" in CLI
count_of_T2=length(grep(x = colnames(f),
                           pattern = "^T2."))
count_of_T2 #--number of samples with Stage 2

#Filtering----
#Filtering genes with an Expression value > 0 in at least two samples out of 44 samples 
filter <- apply(f, 1, function(x) length(x[x>0])>=2)
filtered <- f[filter,]
genes <- rownames(filtered)

t1<-rep(c("T1"),each=count_of_T1)
t2<-rep(c("T2"),each=count_of_T2)
x<-c(t1,t2)
#factor means dividing into levels: T1 and T2
x<-as.factor(x)
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x,
                                                  row.names=colnames(filtered)))
set

#Setting Color theme
colors <- brewer.pal(3, "Set2")

#Plotting basic graphs
#Relative log expression (RLE) plots are a simple, yet powerful, tool for visualizing unwanted variation in high dimensional data
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#Normalize using upper quartile
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


differences <- makeGroups(x)
set3 <- RUVs(set, genes, k=1, differences)
pData(set3)

#Computing Differential expression
dds <- DESeqDataSetFromMatrix(countData = counts(set3),
                              colData = pData(set3),
                              design= ~ W_1 + x)

dds <- DESeq(dds)

res <- results(dds)
head(res)

write.table(res, file = "RUVseq_analysis.csv",
            sep = ",", col.names = NA,
            qmethod = "double")

df <- read.csv("RUVseq_analysis.csv")
#padj: p-val after correction
df <- subset(df, df$padj <= 0.05)

#Volcano plot: differential gene expression
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

#Heatmap
f1 <- f[1:4,]
pheatmap(f1, cluster_cols = FALSE, scale="row")


#Boxplot
f <- t(f)
f <- as.data.frame(f)
View(f[1:10,1:10])

Tumor_stage <- c(rep(c("T1"),each=count_of_T1),
                 rep(c("T2"),each=count_of_T2))
f <- cbind(Tumor_stage,f)
View(f[,1:10])

p<-ggplot(f, aes(x=Tumor_stage, y=C2orf69P3, fill=Tumor_stage)) +
  geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
p

#Violin Plot
p<-ggplot(f, aes(x=Tumor_stage, y=C2orf69P3, fill=Tumor_stage)) +
  geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2))
p
