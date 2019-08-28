


#TO RUN: bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 0.1 > 0.1.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 0.1 > 0.1.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 0.5 > 0.5.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 1 > 1.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 5 > 5.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 10 > 10.Rout'
#bsub -q long -R"select[mem>800] rusage[mem=800]" -M800 -o missensemodel.out 'Rscript --verbose missensemodel_farm_parallel.R 20 > 20.Rout'
#----------FUNCTIONS-------------------

SimMisModel <- function(prophi,emis,df,phi = FALSE,nsim = 10,edist = "gamma",shape = NULL,mixprop = 0.5){
  
  known <- which(df$padj_shetres<0.025 & !is.na(df$padj_shetres))
  unknown <- 1:nrow(df)
  unknown <- unknown[!(unknown %in% known)]
  
  #choose HI genes
  if(phi){
    #prob of being HI function of power to detect
    higenes <- sample(unknown,round(prophi*length(unknown)),prob = df$pHI[unknown])
  } else{
    #all genes have equal prob of being HI
    higenes <- sample(unknown,round(prophi*length(unknown)))
  }
  
  #get enrichment multiplier vector
  mult_mis <- rep(1,nrow(df))
  if(edist == "gamma"){
    mult_mis[higenes] <- 1 + rgamma(length(higenes),scale = emis/shape, shape = shape)
    mult_mis[known] <- df$misratio[known]
  }else if(edist == "mix"){
    nd1 <- round(length(higenes)*mixprop)
    mult_mis[higenes[1:nd1]] <- 1 + rgamma(nd1,scale = emis/shape[1], shape = shape[1])
    mult_mis[higenes[(nd1+1):length(higenes)]] <- 1 + rgamma(length(higenes)-nd1,scale = emis/shape[2], shape = shape[2])
    mult_mis[known] <- df$misratio[known]
  }else{
    mult_mis[higenes] <- emis
    mult_mis[known] <- df$misratio[known]
  }
  mult_mis[mult_mis < 1] <- 1
  #get new vector for probability of LoFs across cohort
  mu_mis <- df$misexpected * mult_mis
  
  #simulate 
  sim = rpois(nsim*length(mu_mis),lambda = mu_mis) 
  sim <- matrix(sim,nrow = nsim,byrow = T)
  #pals <- 1-ppois(t(sim),lambda = unk$misexpected)
  
  total_obs <- sum(df$missense_variant,na.rm = T)
  total_nsig <-  sum(df$missense_variant>df$missig)
  total_misratio <- sum(df$misratio>2)
  
  #probability of observing total DNMs
  total <- rowSums(sim)
  p_total <- dpois(total_obs, lambda = median(total),log = T)
  
  #probability of 1 or more significant genes
  nsig <- colSums(t(sim)>df$missig )
  p_nsig <- dpois(total_nsig,lambda = median(nsig),log = T)
  
  #probability of having x genes with lof enrichment more than bla
  misratio <- t(sim)/df$misexpected
  p_misratio <- dpois(total_misratio, lambda = median(colSums(misratio>2)),log = T)
  
  lik <- p_total + p_nsig + p_misratio
  
  return(lik)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#----------LIBRARY AND PATHS---------------------
#explore corrected
# hard coded paths
mywd = "/nfs/ddd0/jk18/dne/model/missensemodel/"
powpath = "/nfs/ddd0/jk18/dne/power/dne_combined_corrected_bonferroni_2019_05_20_wpower.tab"

library(reshape2)
library(ggplot2)
library(data.table)
library(wesanderson)
#library(MASS)
#library(mclust)
#library(LaplacesDemon)
otherpal <- wes_palette("Darjeeling1")

#script ----------------------------
setwd(mywd)

args = commandArgs(trailingOnly=TRUE)
shape_val = as.numeric(args[1])
print(shape_val)

#read in enrichment results
allhpna <- fread(powpath,sep = "\t",stringsAsFactors = F)
allhpna$pLI.y[is.na(allhpna$pLI.y)] <- median(allhpna$pLI.y,na.rm = T)
allhpna <- allhpna[!is.na(allhpna$misexpected)]
edist <- allhpna$lofratio[!is.na(allhpna$lofratio) & allhpna$padj_shetres<0.025 & !is.na(allhpna$padj_shetres)& allhpna$lofcount>0 & allhpna$consensus_gene]
misdist <- allhpna$misratio[!is.na(allhpna$lofratio) & allhpna$padj_shetres<0.025 & !is.na(allhpna$padj_shetres)& allhpna$consensus_gene &allhpna$me<0.05]
esd <- sd(log(edist))
plimulti <- (sum((allhpna$pLI.y>=0.9 & allhpna$powmed>0.8) & ((!is.na(allhpna$sig) & allhpna$sig )| allhpna$consensus_gene))/sum((allhpna$pLI.y>=0.9 & allhpna$powmed>0.8)))/(sum((allhpna$pLI.y<0.9 & allhpna$powmed>0.8) & ((!is.na(allhpna$sig) & allhpna$sig )| allhpna$consensus_gene))/sum((allhpna$pLI.y<0.9 & allhpna$powmed>0.8)))

#th <- 0.05/nrow(allhpna)
th <- 0.05/18856

unk <- allhpna
unk$pHI <- 1-unk$powmed
unk$pHI <- unk$pHI/sum(unk$pHI)
unk$missig <- qpois(th, lambda = unk$misexpected,lower.tail = F)
unk$misratio[is.na(unk$misratio)] <- 0
unk$missense_variant[is.na(unk$missense_variant)] <- 0


#PARAMETERS


#set variables
obs_total <- sum(unk$lofcount)
prophi_vals <- seq(0,0.35,0.0025)
emis_vals <- seq(0,4.5,0.03)


store_divergence_prob <-array(0,dim = c(length(prophi_vals), length(emis_vals)))
dimnames(store_divergence_prob)[[1]] = round(prophi_vals*nrow(unk))
dimnames(store_divergence_prob)[[2]] = emis_vals


#RUN SIMS
start_time <- Sys.time()
for(i in 1:length(prophi_vals)){
  print(i)
  for(j in 1:length(emis_vals)){
    print(j)  
    store_divergence_prob[i,j] <- mean(replicate(10,SimMisModel(prophi_vals[i],emis_vals[j],unk,phi = T,nsim = 10,edist = "gamma", shape = shape_val)))
  }
}
end_time <- Sys.time()
end_time - start_time


filename = paste(c("Missense_modelresults",args[1],toString(Sys.Date()),"v2.tab"),collapse = "_")
write.table(store_divergence_prob,file = filename,sep = "\t")
filename = paste(c("Missense_modelresults",args[1],toString(Sys.Date()),"v2.rds"),collapse = "_")
saveRDS(store_divergence_prob, file = filename)

