##summarize results with BIC tuning, weighted estimation errors, clustering
library(ggplot2)
library(cowplot)
library(patchwork)
library(pracma)
library(reshape)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(kableExtra)

cluster_results = read.csv("JisstPCA/simulations/ComparativeStudies/NetworkSim/results_cluster.csv", header = TRUE, sep = ",")
cluster_results <- cluster_results[!duplicated(cluster_results),]
cluster_summary  <- cluster_results %>% group_by(q, N, method, factor.name)
cluster_summary <- cluster_summary %>% summarise(RI_mean = mean(RI), RI_sd = sd(RI), ARI_mean = mean(ARI), ARI_sd = sd(ARI), nrep = n())
data.frame(cluster_summary[cluster_summary$q==50 &cluster_summary$N == 20, ])

cluster_summary$method[cluster_summary$method == "iHOSVD"] <- "iHOSVD (oracle)"
cluster_summary$method[cluster_summary$method == "iHOOI"] <- "iHOOI (oracle)"
cluster_subset <- cluster_summary[(cluster_summary$method!="JisstPCA0 (oracle)")&
                                    (cluster_summary$method!="G-JisstPCA0 (oracle)")&
                                    (cluster_summary$method!="G-JisstPCA2 (oracle)")&
                                    (cluster_summary$method!="G-JisstPCA2 (oracle)"),]
cluster_subset$factor.name <- factor(cluster_subset$factor.name, levels = c("sample", "networkX1", "networkX2", "networkY1", "networkY2"))
cluster_subset$method <- factor(cluster_subset$method, levels = 
                                  c("JisstPCA2 (BIC)", "G-JisstPCA2 (BIC)", "JisstPCA0 (BIC)", "G-JisstPCA0 (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)"))

q <- 80
#method_list <- c("JisstPCA2 (BIC)", "G-JisstPCA2 (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)")
method_list <- c("JisstPCA2 (BIC)", "G-JisstPCA2 (BIC)", "JisstPCA0 (BIC)", "G-JisstPCA0 (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)")
for(N in c(20, 40)){
  results0 <- NULL
  for(method in method_list){
    ind <- which(cluster_subset$q==q & cluster_subset$N==N & cluster_subset$method == method)
    results_tmp <- cluster_subset[ind, ]
    results_tmp <- results_tmp[order(results_tmp$factor.name),c("ARI_mean", "ARI_sd")]
    results_tmp <- round(results_tmp, 3)
    results_tmp <- apply(results_tmp, 1, function(x){paste(x[1], "(", x[2], ")", sep = "")}) 
    results0 <- cbind(results0, results_tmp)
  }
  if(N == 20){
    table <- results0
  }else{
    table <- cbind(table, results0)
  }
}
table <- data.frame(table)
row.names(table) <- NULL; colnames(table) <- NULL
kable(table, format = "latex", booktabs = T) %>%
  add_header_above(c("N=20" = 4, "N=40" = 4)) 

results = read.csv("JisstPCA/simulations/ComparativeStudies/NetworkSim/results.csv", header = TRUE, sep = ",")
results_summary  <- results %>% group_by(q, N, method, factor.name, factor.number)
results_summary <- results_summary %>% summarise(err_mean = mean(err), err_sd = sd(err), nrep = n())
results_summary <- results_summary %>% unite("factor", 
                                             c(factor.name,factor.number), 
                                             remove = FALSE)
results_summary$factor <- factor(results_summary$factor, levels = c("u_1", "u_2", "V_1", "V_2", "W_1", "W_2"))
results_summary$method[results_summary$method == "iHOSVD"] <- "iHOSVD (oracle)"
results_summary$method[results_summary$method == "iHOOI"] <- "iHOOI (oracle)"

#results_summary$method <- factor(results_summary$method, levels = c("JisstPCA", "G-JisstPCA", "iHOSVD", "iHOOI"))
#method_list <- c("JisstPCA2 (BIC)", "G-JisstPCA2 (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)")
method_list <- c("JisstPCA2 (BIC)", "G-JisstPCA2 (BIC)", "JisstPCA0 (BIC)", "G-JisstPCA0 (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)")
results_sub <- results_summary[results_summary$method %in% method_list, ]

q <- 80
  for(N in c(20, 40)){
    for(method in method_list){
      table_tmp <- results_sub[results_sub$q == q & results_sub$N == N & results_sub$method == method, 
                               c("factor", "err_mean", "err_sd")]
      table_tmp <- table_tmp[order(table_tmp$factor), c(2, 3)]
      table_tmp <- round(table_tmp, 3)
      table_tmp[,3] <- apply(table_tmp, 1, function(x){paste(x[1], "(", x[2], ")", sep = "")}) 
      if(method != "JisstPCA2 (BIC)"){
        table0 <- cbind(table0, table_tmp[,3])
      }else{
        table0 <- table_tmp[,3]
      }
    }
    table0 <- rbind(method_list, table0)
    if(N == 20){
      table <- table0
    }else{
      table <- cbind(table, table0)
    }
  }

row.names(table) <- NULL; colnames(table) <- NULL
kable(table, format = "latex", booktabs = T) %>%
  add_header_above(c("N=20" = 3, "N=40" = 3)) 


