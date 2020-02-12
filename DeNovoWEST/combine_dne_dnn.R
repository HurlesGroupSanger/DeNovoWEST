# Joanna Kaplanis and Kaitlin Samocha
#
#R script that takes output from enrichment test DeNovoWEST on all variants, just missense variants, and a missense clustering test from DeNovoNear and combines p-values
#
# Usage: Rscript combine_IHWcorrect.R [dne_all_file] [dne_mis_file] [dnn_file] [out_file]
#
#
#-----------------IMPORTS--------------
library(metap)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
dne_all_file = args[1]
dne_mis_file = args[2]
dnn_file = args[3]
out_file = args[4]

#dnn_file = "/Volumes/ddd0/jk18/dne/DeNovoWEST/input/denovonear_out_missense_31058_ntrios_2019_05_15.txt"
#dne_file_all = "/Volumes/ddd0/jk18/dne/newtest/dne_test_2020_01_21_all.tab"
#dne_file_mis = "/Volumes/ddd0/jk18/dne/newtest/missense_test/dne_test_2020_01_21_mis.tab"
#out_file = "/Volumes/ddd0/jk18/dne/newtest/dne_combined_corrected_bonferroni_2020_02_12.tab"

#dnn_file = "/Volumes/ddd0/ks20/active_metaDNM_analyses/caddv1.0_4sets/undiagnosed_analyses/denovonear_out_missense_undiagnosed_24288_ntrios_2020_02_11.txt"
#dne_file_all = "/Volumes/ddd0/jk18/dne/newtest/undiagnosed/dne_test_2020_02_10_ud.tab"
#dne_file_mis = "/Volumes/ddd0/jk18/dne/newtest/undiagnosed/missense/dne_test_mis_2020_02_10_ud.tab"
#out_file = "/Volumes/ddd0/jk18/dne/newtest/dne_combined_corrected_bonferroni_2020_02_12_ud.tab"


#------COMBINE RESULTS----------------------

#read in de novo near output
dnn <- fread(dnn_file,sep = "\t",header = T, stringsAsFactors = F)
dnn$dnn_p <- dnn$probability

#read in enrichment for all variants output
dne_all <- fread(dne_file_all,sep = "\t",header = T,stringsAsFactors = F)
names(dne_all)  <- c("symbol","hgnc_id","expected_all","observed_all","dne_all_p","info_all")

dne_mis <- fread(dne_file_mis,sep = "\t",header = T,stringsAsFactors = F)
names(dne_mis) <- c("symbol","hgnc_id","expected_mis","observed_mis","dne_mis_p","info_mis")

#merge data
res <-merge(dne_all,dnn,by.x = "symbol",by.y = "gene_id",all = T)
res<- merge(res,dne_mis,by = c("symbol","hgnc_id"),all = T)

# sumlog combines p values
res$com_p <- NA
res$com_p[(!is.na(res$dne_mis_p) & !is.na(res$dnn_p))] <- apply(res[(!is.na(res$dne_mis_p) & !is.na(res$dnn_p)),c("dnn_p","dne_mis_p")],MARGIN = 1, FUN = function(pvals) sumlog(pvals)$p)

# take min of combined and of enrichment
res$min_p <- apply(res[,c("com_p","dne_all_p")],MARGIN = 1, FUN = min,na.rm = T)
res$min_p[which(!is.finite(res$min_p))] <- 1

#write to file
write.table(res,file = out_file,sep = "\t",col.names = T,row.names = F,quote = F)


