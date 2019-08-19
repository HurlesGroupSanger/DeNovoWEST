# Joanna Kaplanis and Kaitlin Samocha
#
#R script that takes output from enrichment test DeNovoWEST and clustering test from DeNovoNear and combines p-values
# and applies IHW weighting based on s_het corrected for mutability
#
# Usage: Rscript combine_IHWcorrect.R [dne_file] [dnn_file] [out_file]
#
#
#--------------------------------------

# filepaths and set up
mut_file = 'input/gene_mutability_length.tab'
hgnc_genefile = "input/protein-coding_gene.txt"
shet_file = 'input/s_het.csv'

#-----------------IMPORTS--------------
source("useful_functions.R")
mypal<-wes_palette("Royal1")
library(lattice)
library(metap)
library(ggplot2)
library(data.table)
library(gridExtra)
library(wesanderson)
#devtools::install_github("nignatiadis/IHW")
library(IHW)

args = commandArgs(trailingOnly=TRUE)
dne_file = args[1]
dnn_file = args[2]
out_file = args[3]

#---------LOAD DATA--------------------
# first part is to load the HGNC genes so that all other gene sets can be annotated with them
hgnc_genes <- fread(hgnc_genefile,sep = "\t",header = T, stringsAsFactors = F)
# to get rid of the HGNC: part of the hgnc_id column
hgnc_genes$num_hgnc_id <- apply(hgnc_genes, 1, function(row) strsplit(row[1], ':')[[1]][2])
hgnc_genes$chrh <-  sapply(strsplit(sapply(strsplit(hgnc_genes$location,split = "q"), "[", 1),"p"),"[",1)

#mutability info
muts <- read.table(mut_file,header = T, stringsAsFactors = F, sep = "\t")
muts$glen <- muts$pos

#S_HET
shet <- read.csv(shet_file,sep = ",",header = T,stringsAsFactors = F)
# add a column that lists NA if the gene_symbol isn't in HGNC
shet$matchedhgnc <- apply(shet, 1, function(row) if(row[1] %in% hgnc_genes$symbol) { row[1] } else { NA } )
shet$newsymbol <- shet$gene_symbol
shet$newsymbol[shet$gene_symbol%in%names(myl)] <- myl[shet$gene_symbol[shet$gene_symbol%in%names(myl)]]

# use matchedhgnc (row[16])
# if NA, use newsymbol (row[17])
shet$replacementsymbol <- apply(shet, 1, function(row) if(is.na(row[16])){ row[17] } else { row[16] })

# clean up
shet$symbol <- shet$replacementsymbol
shet$matchedhgnc <- NULL
shet$newsymbol <- NULL
shet$replacementsymbol <- NULL

#---------MERGE--------------------------
genes <- merge(muts,shet, by.x = "symbol",by.y = "symbol",all.x = T)
genes <- merge(genes,hgnc_genes,by.x = "symbol",by.y = "symbol",all.y = T)
genes$chr <- genes$chrh
genes <- genes[!is.na(genes$prob),]

#---------GET RESIDUALS--------------------

mshetX <- median(genes$s_het[genes$chr == "X" & !is.na(genes$chr)],na.rm = T)
msheta <- median(genes$s_het[genes$chr != "X" & !is.na(genes$chr)],na.rm = T)
genes$s_het[is.na(genes$s_het) & (genes$chr != "X" | is.na(genes$chr))] <- msheta
genes$s_het[is.na(genes$s_het) & genes$chr == "X"] <- mshetX

#get shet residuals after regressing out gene length and mutability
shetgl <- lm(s_het~prob,data = genes)
genes$shetres <- shetgl$residuals

#------COMBINE RESULTS----------------------

#read in de novo near output
dnn <- fread(dnn_file,sep = "\t",header = T, stringsAsFactors = F)
dnn$dnn_p <- dnn$probability

#read in enrichment output
dne <- fread(dne_file,sep = "\t",header = T,stringsAsFactors = F)
dne$dne_p <- dne$`p-value`

#merge data
res <-merge(dne,dnn,by.x = "symbol",by.y = "gene_id",all = T)
genes <- merge(genes,res,by.x = "symbol",by.y = "symbol",all = T)
genes <- genes[!is.na(genes$prob),]

# sumlog combines p values
genes$com_p <- apply(genes[,c("dnn_p","dne_p")],MARGIN = 1, FUN = min,na.rm = T)
genes$com_p[(!is.na(genes$dne_p) & !is.na(genes$dnn_p))] <- apply(genes[(!is.na(genes$dne_p) & !is.na(genes$dnn_p)),c("dnn_p","dne_p")],MARGIN = 1, FUN = function(pvals) sumlog(pvals)$p)

# take min of combined and of enrichment
genes$min_p <- apply(genes[,c("com_p","dne_p")],MARGIN = 1, FUN = min,na.rm = T)
genes$min_p[which(!is.finite(genes$min_p))] <- 1
genes <- genes[!is.na(genes$chr),]

#-------APPLY IHW-------------------

#set threshold
mya <- 0.025
#correction type
atype <- "bonferroni"
#number of genes
ngenes <- length(genes$symbol)
#IHW bins
nbins <- 7

#get IHW bins for s_het residuals
genes$shetresb <- groups_by_filter(genes$shetres,nbins)
m_groups_shetres <- table(genes$shetresb)

#apply IHW
ihwshetres <- ihw(min_p ~ shetresb,  data = tres, alpha = mya,adjustment_type = atype,m_groups = m_groups_shetres)
#save p-values
tres$padj_shetres <- (adj_pvalues(ihwshetres))

#write to file
write.table(tres, out_file, quote = F, row.names = F, sep = "\t")
