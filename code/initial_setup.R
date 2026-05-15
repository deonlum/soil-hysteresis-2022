## Set-up source file for analysis

## R script for: 
rm(list = ls())

# 1. Set-up ====
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggnewscale)
library(gridExtra)
library(reshape2)
library(cowplot)
library(mgcv)
library(metR) # for contour text labels

mycols = c("field" = "#316294", 
           "drought" = "#EC2124")
theme_set(theme_bw())

## 1.1 Reading in main dataset ====
main_df = read.csv("./data/Dataset_S1.csv", stringsAsFactors = TRUE)

## Ordering factors
# timepoints
main_df$timepoint = factor(main_df$timepoint, 
                           levels = c("3 days", "7 days", "14 days",
                                      "35 days", "70 days", "c1", "c14"))

# treatments
main_df$treatment = factor(main_df$treatment,
                           levels = c("control", "control-field",
                                      "control-drought", "field",
                                      "drought"))

# Calculating Soil Water Deficit (SWD) as the fraction of soil water content
# relative to the water the soil can hold 
# (i.e., its water holding capacity, measured as 54.6605%)
main_df$swd = (1-main_df$wc/54.6605)

## Assigning the order for joining points later
main_df$point_order = c(c(12:22, 1:11),
                        c(12:22, 1:11)+22,
                        c(12:22, 1:11)+44,
                        c(12:22, 1:11)+66,
                        c(12:22, 1:11)+88,
                        c(0,0,0,0,0,0,0,0,0))

## 1.2 Loading in sequencing datasets ====
bac_seq_data = read.csv("./data/Dataset_S2.csv")
fun_seq_data = read.csv("./data/Dataset_S3.csv")

dim(bac_seq_data) # sequences are in the first column, with taxonomy from 121:127

## Getting taxonomy
bac_taxonomy = bac_seq_data[,c(121:127)]
fun_taxonomy = fun_seq_data[,c(121:127)]
rownames(bac_taxonomy) = paste0("ASV", 1:nrow(bac_taxonomy))
rownames(fun_taxonomy) = paste0("ASV", 1:nrow(fun_taxonomy))

## Getting counts, transpose for use with vegan
bac_counts = t(as.matrix(bac_seq_data[,c(2:120)]))
fun_counts = t(as.matrix(fun_seq_data[,c(2:120)]))
rownames(bac_counts) = 1:119
colnames(bac_counts) = paste0("ASV", 1:ncol(bac_counts))
rownames(fun_counts) = 1:119
colnames(fun_counts) = paste0("ASV", 1:ncol(fun_counts))

## 1.3 Removing low abundance reads from sequencing data ====
bac_lib_sizes = apply(bac_counts,1,sum)
fun_lib_sizes = apply(fun_counts,1,sum)

# Removing samples that are never higher than 0.01% of a sample
# This roughly translates to fewer than ~3-5 reads
bac_asv_sizes = apply(bac_counts, 2, sum)
fun_asv_sizes = apply(fun_counts, 2, sum)

remove_rare_ASV = function(my_matrix){
  lib_size = apply(my_matrix, 1, sum) # Find the library size
  low_prop = apply(my_matrix, 2, function(x){
    sum(x>(lib_size/10000)) # looking for ASVs that are in less than 0.01% of a sample
  })
  # Remove ASV if they never make up more than 0.01% of any sample
  new_matrix = my_matrix[,low_prop != 0] 
}
clean_bac = remove_rare_ASV(bac_counts)
clean_fun = remove_rare_ASV(fun_counts)

## 1.4 Removing taxa of non-interest ====
bac_tokeep = rownames(bac_taxonomy)[which(bac_taxonomy$Kingdom == "Bacteria" |
                                            bac_taxonomy$Kingdom == "Archaea")]
fun_tokeep = rownames(bac_taxonomy)[which(fun_taxonomy$Kingdom == "k__Fungi")]
clean_bac = clean_bac[,colnames(clean_bac) %in% bac_tokeep]
clean_fun = clean_fun[,colnames(clean_fun) %in% fun_tokeep]

# Calculating new library sizes after removing ASVs
clean_bac_lib = apply(clean_bac, 1, sum)
clean_fun_lib = apply(clean_fun, 1, sum)

# Making sure taxonomy reflects new dataset
clean_bac_taxonomy = bac_taxonomy[rownames(bac_taxonomy) %in% colnames(clean_bac),]
clean_fun_taxonomy = fun_taxonomy[rownames(fun_taxonomy) %in% colnames(clean_fun),]


## Functions
get_gam_ci = function(my_mod, new_data){
  
  # Adapted from: https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
  rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
  }
  
  # Pointwise CI
  Vb <- vcov(my_mod)
  newd <- new_data
  pred <- predict(my_mod, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  
  # Simulated CI
  N <- 10000
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  Cg <- predict(my_mod, newd, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  return(data.frame(fit = pred$fit,
                    point_lo = pred$fit - 2*se.fit,
                    point_hi = pred$fit + 2*se.fit,
                    sim_lo = pred$fit - crit*se.fit,
                    sim_hi = pred$fit + crit*se.fit))
  
}

## Repeated rarefaction for diversity indices
repeat_rarefy_div = function(my_seq, nrep = 100){
  
  # Finding smallest library
  lib_sizes = apply(my_seq, 1, sum)
  
  # Creating objects to store values
  avg_richness = numeric(nrow(my_seq))
  avg_shannon = numeric(nrow(my_seq))
  avg_matrix = matrix(0, nrow = nrow(my_seq), ncol = nrow(my_seq))
  
  # Repeatedly rarefying
  for (i in 1:nrep){
    if (i%%5 == 0){
      cat("running iteration", i, "\n")
    }
    
    ## Rarefy
    temp_rare = rrarefy(my_seq, sample = min(lib_sizes))
    
    ## Getting diversity values
    temp_richness = apply(temp_rare, 1, function(x) sum(x>0))
    temp_shannon = vegan::diversity(temp_rare, index = "shannon")
    temp_dist = vegdist(decostand(temp_rare, method = "hellinger"),
                        method = "bray")
    
    # Calculating averages and storing values
    avg_richness = avg_richness + temp_richness/nrep
    avg_shannon = avg_shannon + temp_shannon/nrep
    avg_matrix = avg_matrix + as.matrix(temp_dist)/nrep
  }
  return(list(avg_matrix = avg_matrix, 
              avg_richness = avg_richness, 
              avg_shannon = avg_shannon))
}


## Repeated rarefaction of reads
repeat_rarefy = function(my_seq, nrep = 100){
  
  # Finding smallest library
  lib_sizes = apply(my_seq, 1, sum)
  
  # Creating objects to store values
  avg_matrix = matrix(0, nrow = nrow(my_seq), ncol = ncol(my_seq))
  
  # Repeatedly rarefying
  for (i in 1:nrep){
    if (i%%5 == 0){
      cat("running iteration", i, "\n")
    }
    
    ## Rarefy
    temp_rare = rrarefy(my_seq, sample = min(lib_sizes))
    
    ## Getting diversity values
    avg_matrix = avg_matrix + temp_rare/nrep
  }
  return(avg_matrix)
}
