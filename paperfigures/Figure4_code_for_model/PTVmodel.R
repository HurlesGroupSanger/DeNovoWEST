#
#Author: Joanna Kaplanis @queenjobo 
#Parts based on code from Patrick Short @pjshort
#16/08/2019
#
# Usage: Rscript --verbose PTVmodel.R > PTVmodel.Rout'
#

# Simulation Function
SimModel <- function(prophi,elof,df,phi = FALSE,nsim = 10,edist = FALSE,esd = NULL){
  
  known <- which(df$sig | df$consensus_gene)
  resnown <- 1:nrow(df)
  resnown <- resnown[!(resnown %in% known)]
  
  #choose HI genes
  if(phi){
    #prob of being HI function of power to detect
    higenes <- sample(resnown,round(prophi*length(resnown)),prob = res$pHI[resnown])
  } else{
    #all genes have equal prob of being HI
    higenes <- sample(resnown,round(prophi*length(resnown)))
  }
  
  #get enrichment multiplier vector
  mult_lof <- rep(1,nrow(res))
  if(edist){
    mult_lof[higenes] <- exp(rnorm(length(higenes),mean = log(elof), sd = esd))
    mult_lof[known] <- df$lofratio[known]
  }else{
    mult_lof[higenes] <- elof
    mult_lof[known] <- df$lofratio[known]
  }
  mult_lof[mult_lof<1] <- 1
  #get new vector for probability of LoFs across cohort
  mu_lof <- df$lofexpected * mult_lof
  
  #simulate 
  sim = rpois(nsim*length(mu_lof),lambda = mu_lof) 
  sim <- matrix(sim,nrow = nsim,byrow = T)
  #pals <- 1-ppois(t(sim),lambda = res$lofexpected)
  
  total_obs <- sum(df$lofcount)
  total_nsig <- sum(df$lofcount>df$lofsig)
  total_lofratio <- sum(df$lofratio>2)
  
  #probability of observing total DNMs
  total <- rowSums(sim)
  p_total <- dpois(total_obs, lambda = median(total),log = T)
  
  #probability of 1 or more significant genes
  nsig <- colSums(t(sim)>df$lofsig )
  p_nsig <- dpois(total_nsig,lambda = median(nsig),log = T)
  
  #probability of having x genes with lof enrichment more than bla
  lofratio <- t(sim)/df$lofexpected
  p_lofratio <- dpois(total_lofratio, lambda = median(colSums(lofratio>2)),log = T)
  
  lik <- p_total + p_nsig + p_lofratio
  
  return(lik)
}

#------PATHS----------------
powpath = "../../input/extended_denovoWEST_results.tab"
#import libraries-------------------
library(reshape2)
library(ggplot2)
library(data.table)
library(wesanderson)

otherpal <- wes_palette("Darjeeling1")

#script ----------------------------
set.seed(5)

#read in enrichment results
res <- fread(powpath,sep = "\t",stringsAsFactors = F)
#replace genes without pLI with median
res$pLI[is.na(res$pLI)] <- median(res$pLI,na.rm = T)
res <- res[!is.na(res$lofexpected)]

edist <- res$lofratio[!is.na(res$lofratio) & res$sig & res$lofcount>0 & res$consensus_gene]
esd <- sd(log(edist))

#For well powered genes, how much more likely are you to be a DD gene if you have pLI>0.9, this is fed into the model
plimulti <- (sum((res$pLI>=0.9 & res$powmed>0.8) & ((!is.na(res$sig) & res$sig )| res$consensus_gene))/sum((res$pLI>=0.9 & res$powmed>0.8)))/(sum((res$pLI<0.9 & res$powmed>0.8) & ((!is.na(res$sig) & res$sig )| res$consensus_gene))/sum((res$pLI<0.9 & res$powmed>0.8)))

#significance threshold 
ngenes <- 18762
th <- 0.05/ngenes

# set the probability of being haploinsufficient as a function of power and pLI
res$pHI <- 1-res$powmed
res$pHI[res$pLI>=0.9] <- (1-res$powmed[res$pLI>=0.9])*plimulti
res$pHI <- res$pHI/sum(res$pHI)

#calculate the number of PTVs needed to cross significance threshold for each gene
res$lofsig <- qpois(th, lambda = res$lofexpected,lower.tail = F)
#the total number of observed PTVs
obs_total <- sum(res$lofcount)

#set variables
#range of proportion of genes that are haploinsufficient DD genes
prophi_vals <- seq(0.0001,0.20,0.0025)
#range of PTV enrichment values
elof_vals <- seq(0,25,0.25)

#initialise matrix
store_divergence_prob <-matrix(ncol=length(prophi_vals), nrow=length(elof_vals))
colnames(store_divergence_prob) = round(prophi_vals*nrow(res))
rownames(store_divergence_prob) = elof_vals

#conduct simulation across varying values
start_time <- Sys.time()
for(i in 1:length(prophi_vals)){
  for(j in 1:length(elof_vals)){
    store_divergence_prob[j,i] <- mean(replicate(30,SimModel(prophi_vals[i],elof_vals[j],res,phi = T,nsim = 20,edist = F, esd = esd)))
  }
}
end_time <- Sys.time()
end_time - start_time

m = melt(store_divergence_prob, varnames = c("elof_vals", "prophi_vals"), value.name = "logprob")
m$prob <- exp(m$logprob)

filename = paste(c("PTV_modelresults_",toString(Sys.Date()),".tab"),collapse = "")
write.table(m,file = filename,sep = "\t",row.names = F,col.names = T,quote = F)
