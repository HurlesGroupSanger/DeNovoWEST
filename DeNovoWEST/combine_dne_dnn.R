# Joanna Kaplanis and Kaitlin Samocha
#
#R script that takes output from enrichment test DeNovoWEST on all variants, just missense variants, and a missense clustering test from DeNovoNear and combines p-values
#
# Usage: Rscript combine_dne_dnn.R [dne_all_file] [dne_mis_file] [dnn_file] [out_file]
#
#
#-----------------IMPORTS--------------
library(metap)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
dne_file_all= args[1]
dne_file_mis = args[2]
dnn_file = args[3]
out_file = args[4]

#dnn_file = "../input/denovonear_out_missense_31058_ntrios_2019_05_15.txt"
#dne_file_all = "../input/merged_all_dne_test_ppv_2020_03_09.tab"
#dne_file_mis = "../input/merged_mis_dne_test_ppv_2020_03_09.tab"
#out_file = "../input/dne_combined_2020_03_09.tab"
#------COMBINE RESULTS----------------------

#read in de novo near output
dnn <- fread(dnn_file,sep = "\t",header = T, stringsAsFactors = F)
dnn$dnn_p <- dnn$probability
dnn <- dnn[dnn$mutation_category == "missense"]

#read in enrichment for all variants output
dne_all <- fread(dne_file_all,sep = "\t",header = T,stringsAsFactors = F)
names(dne_all)  <- c("symbol","hgnc_id","expected_all","observed_all","dne_all_p","info_all")

dne_mis <- fread(dne_file_mis,sep = "\t",header = T,stringsAsFactors = F)
names(dne_mis) <- c("symbol","hgnc_id","expected_mis","observed_mis","dne_mis_p","info_mis")

#merge data
res <-merge(dne_all,dnn,by.x = "symbol",by.y = "gene_id",all.x = T)
res<- merge(res,dne_mis,by = c("symbol","hgnc_id"),all = T)

# sumlog combines p values
res$com_p <- NA
res$com_p[(!is.na(res$dne_mis_p) & !is.na(res$dnn_p))] <- apply(res[(!is.na(res$dne_mis_p) & !is.na(res$dnn_p)),c("dnn_p","dne_mis_p")],MARGIN = 1, FUN = function(pvals) sumlog(pvals)$p)

# take min of combined and of enrichment
res$min_p <- apply(res[,c("com_p","dne_all_p")],MARGIN = 1, FUN = min,na.rm = T)
res$min_p[which(!is.finite(res$min_p))] <- 1

#write to file
write.table(res,file = out_file,sep = "\t",col.names = T,row.names = F,quote = F)

