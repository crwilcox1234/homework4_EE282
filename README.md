Homework 4 answers
======================

Part1 : Summarize partitions of a genome assembly
-----------------------------------------------------

Answers and code for this section will be in the Partitions_genome_assembly directory
-----------------------------------------------------------------------------------------

For this section the Drosophila melanogaster genome was used: dmel-all-chromosome-r6.24.fasta

bioawk was used to partition the genome into over 100kb and under 100kb:

```bash

bioawk -c fastx 'length($seq) > 100000{ print ">"$name; print $seq }'  dmel-all-chromosome-r6.24.fasta > dmel-all-chr-over100kb.fasta

bioawk -c fastx 'length($seq) <= 100000{ print ">"$name; print $seq }'  dmel-all-chromosome-r6.24.fasta > dmel-all-chr-less_eq_100kb.fasta

```

##Calculate the following for Whole genome

```bash

faSize dmel-all-chromosome-r6.24.fasta

```

1. Total number of nucleotides: 143726002 bases
2. Total number of Ns: 1152978 N's
3. Total number of sequences: 1870 sequences

##Calculate the following for >100kb:

```bash

faSize dmel-all-chr-over100kb.fasta

```

1. Total number of nucleotides: 137547960  bases
2. Total number of Ns: 490385  N's
3. Total number of sequences: 7 sequences


##Calculate the following for <= 100kb:

```bash

faSize dmel-all-chr-less_eq_100kb.fasta

```

1. Total number of nucleotides: 6178042  bases
2. Total number of Ns: 662593  N's
3. Total number of sequences: 1863 sequences


Plots of the following for the whole genome, for all sequences <= 100kb, and all sequences > 100kb:
---------------------------------------------------------------------------------------------------------

1. Sequence length distribution

###Will use this file to plot both sequence length distribution and cumulative genome size sorted from largest to smallest

The same code was used for whole genome, >100kb, and <=100kb, changing the input file to the files containing reads from 
the whole genome, >100kb, and <=100kb

```bash

bioawk -c fastx '{ print $name, length($seq) }' dmel-all-chromosome-r6.24.fasta > Sequence_length.txt

```

2. Sequence GC% distribution

##The GC content:

The same code was used for whole genome, >100kb, and <=100kb, changing the input file to the files containing reads from
the whole genome, >100kb, and <=100kb

```bash

bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta > GC_content.txt

```

3.  Cumulative genome size sorted from largest to smallest sequences

The sequence length distribution file for each genome size file (whole, >100kb, and <=100kb) will be used to plot both sequence length 
distribution and cumulative genome size sorted from largest to smallest

Alt-PLOTS
=========

R code for plots. This code includes the code for the 9 plots. The plots for each question have the same scale, except for 3 additional plots
(for a total of 12 plots) in the >100kb section.  Since there are only 7 sequences over 100kb, the scale of the extra 3 plots was reduced.

```

#install.packages("ggplot2")
library(ggplot2)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

##Whole GENOME
setwd("D:/DATA/UC Irvine/JJEMERSON_CLASS_2018/homework4/whole")
dmel <- read.table("GC_content.txt", header = FALSE)

#qplot(dmel_l$V2, geom="histogram") 


a <- ggplot(data=dmel, aes(dmel$V2)) + geom_histogram(
  bins = 10 ,
  fill="red",alpha = .5) +
  labs(title="Whole_Genome_GC_Content") +
  labs(x="GC_Content", y="Count") 

print(a)


dmel_l <- read.table("Sequence_length.txt", header = F, row.names=1)
log <- log10(dmel_l)

b <- ggplot(data=log, aes(log$V2)) + geom_histogram(
  bins = 10 ,
  #binwidth = 1000
  fill="green", alpha=0.5) +
  labs(title="Length_Distribution_Whole") +
  labs(x="sequence_length_base(log10)", y="Count") 

print(b)

#sort in decending order
data <- as.data.frame(dmel_l[order(-dmel_l$V2),])
colnames(data)[1] <- "V2"
#data <- dmel_l[order(-dmel_l$V2),]
#Cumulative genome size 
data$n_seqs <- 1:nrow(data)
data$gen_size <- cumsum(data$V2)

c <- ggplot(data, aes(n_seqs, gen_size)) + geom_point()+ geom_line()

print(c)

##GENOME over 100kb
setwd("D:/DATA/UC Irvine/JJEMERSON_CLASS_2018/homework4/over100kb")
dmel_over <- read.table("GC_content.txt", header = FALSE)
d <- ggplot(data=dmel_over, aes(dmel_over$V2)) + geom_histogram(
  bins = 10 ,
  fill="purple",alpha = .5) +
  labs(title="GC_Content_over100kb") +
  labs(x="GC_Content", y="Count") +
  coord_cartesian(xlim=c(0,0.75),ylim=c(0,700))

print(d)

e <- ggplot(data=dmel_over, aes(dmel_over$V2)) + geom_histogram(
  bins = 10 ,
  fill="purple",alpha = .5) +
  labs(title="GC_Content_over100kb") +
  labs(x="GC_Content", y="Count") +
  coord_cartesian(xlim=c(0,0.75),ylim=c(0,5))

print(e)

dmel_over_l <- read.table("Sequence_length.txt", header = F, row.names=1)
log_over <- log10(dmel_over_l)

f <- ggplot(data=log_over, aes(log_over$V2)) + geom_histogram(
  bins = 10 ,
  #binwidth = 1000
  fill="orange", alpha=0.5) +
  labs(title="Length_Distribution_over100kb") +
  labs(x="sequence_length_base(log10)", y="Count")+
  coord_cartesian(xlim=c(0,9),ylim=c(0,1500))

print(f)

g <- ggplot(data=log_over, aes(log_over$V2)) + geom_histogram(
  bins = 10 ,
  #binwidth = 1000
  fill="orange", alpha=0.5) +
  labs(title="Length_Distribution_over100kb") +
  labs(x="sequence_length_base(log10)", y="Count")+
  coord_cartesian(xlim=c(0,9),ylim=c(0,7))

print(g)

##sort in decending order
data_over_l <- dmel_over_l[order(-dmel_l$V2),]
#Cumulative genome size 
data_over_l <- as.data.frame(dmel_over_l[order(-dmel_over_l$V2),])
colnames(data_over_l)[1] <- "V2"
#data <- dmel_l[order(-dmel_l$V2),]
#Cumulative genome size 
data_over_l$n_seqs <- 1:nrow(data_over_l)
data_over_l$gen_size <- cumsum(data_over_l$V2)

h <- ggplot(data_over_l, aes(n_seqs, gen_size)) + geom_point()+ geom_line()+
  coord_cartesian(xlim=c(0,2000),ylim=c(0,150000000))

print(h)

##GENOME under 100kb
setwd("D:/DATA/UC Irvine/JJEMERSON_CLASS_2018/homework4/under100kb")
dmel_under <- read.table("GC_content.txt", header = FALSE)
i <- ggplot(data=dmel_under, aes(dmel_under$V2)) + geom_histogram(
  bins = 10 ,
  fill="blue",alpha = .5) +
  labs(title="GC_content_under100kb") +
  labs(x="GC_Content", y="Count")
  #coord_cartesian(xlim=c(0,0.75),ylim=c(0,700))

print(i)

dmel_under_l <- read.table("Sequence_length.txt", header = F, row.names=1)
log_under <- log10(dmel_under_l)

k <- ggplot(data=log_under, aes(log_under$V2)) + geom_histogram(
  bins = 4 ,
  #binwidth = 1000
  fill="yellow", alpha=0.5) +
  labs(title="Length_Distribution_under100kb") +
  labs(x="sequence_length_base(log10)", y="Count")+
  coord_cartesian(xlim=c(0,9),ylim=c(0,1500))

print(k)

#sort in decending order
data_under_l <- dmel_under_l[order(-dmel_under_l$V2),]
#Cumulative genome size 
data_under_l <- as.data.frame(dmel_under_l[order(-dmel_under_l$V2),])
colnames(data_under_l)[1] <- "V2"
#data <- dmel_l[order(-dmel_l$V2),]
#Cumulative genome size 
data_under_l$n_seqs <- 1:nrow(data_under_l)
data_under_l$gen_size <- cumsum(data_under_l$V2)

l <- ggplot(data_under_l, aes(n_seqs, gen_size)) + geom_point()+ geom_line()+
  coord_cartesian(xlim=c(0,2000),ylim=c(0,150000000))

print(l)

m <- ggarrange(a, b, c, d, e, f, g, h, i, k, l,
          ncol = 2, nrow = 2)

ggexport(m)

```

Plots in PDF format
-----------------------

To make the plots ggplot2 was used using the geom_histogram function, the geom_point, and geom_line functions in ggplot2.  ggarrange and ggexport were used to arrange all of the plots and export to a single pdf document. The code is in the .R file Plots.R

All plots are in one pdf file: Plots_over_under_whole.pdf

Genome assembly
===================


# homework4_EE282
