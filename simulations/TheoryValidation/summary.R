library(ggplot2)
library(latex2exp)
library(dplyr)
library(grid)
library(gridExtra)
setwd("JisstPCA/simulations/TheoryValidation")
results_file <- "JisstPCA/simulations/TheoryValidation/JisstPCA_res.csv"
#results_file <- "JisstPCA/simulations/TheoryValidation/G-JisstPCA_res.csv"
fig_path <- "JisstPCA/simulations/TheoryValidation/fig/Jisst/"
#fig_path <- "/Applications/Files/research/JisstPCA/JisstPCA_new/validation/fig/GJisst/"

## validation of JisstPCA
results = read.csv(results_file, header = TRUE)
results_iter_summary  <- results %>% group_by(factor, dim, SNR, iteration)
stat_err_results <- results_iter_summary %>% summarise(est_err_mean = mean(stat_err), est_err_sd = sd(stat_err))

results_summary  <- results %>% group_by(dim, SNR, iteration)
loss_results <- results_summary %>% summarise(loss_mean = mean(log(optim_loss)), loss_sd = sd(log(optim_loss)))
loss_results <- loss_results[loss_results$iteration != 0, ]
#optimization loss vs. iteration
for(dim in unique(loss_results$dim)){
  for(SNR in unique(loss_results$SNR)){
    g <- ggplot(loss_results[loss_results$dim==dim&loss_results$SNR == SNR,]) + geom_line(aes(x = iteration, y = loss_mean)) +
      ylab("log(integrated squared loss)") + xlab("Iteration") + scale_x_continuous(breaks = 1:10)
    ggsave(g, file = paste(fig_path, sprintf("optimization_converg_dim%d_SNR%d_2.png", dim, SNR/0.375), sep = ""), 
           width = 8, height = 10, units = "cm")
  }
}
## when SNR is large, optimization loss even increases a little bit (spectral initialization is the best, 
## while further iterations make it slightly worse)

#statistical error vs. iteration]
for(SNR in unique(loss_results$SNR)){
  g <- ggplot(stat_err_results[stat_err_results$SNR == SNR,], aes(x = iteration, y = est_err_mean, color = factor)) + 
    geom_errorbar(aes(ymin = est_err_mean - est_err_sd / sqrt(20) * qnorm(0.975), ymax = est_err_mean + est_err_sd / sqrt(20) * qnorm(0.975)), width=.2) +
    geom_line(aes(x = iteration, y = est_err_mean, color = factor, linetype = as.factor(dim))) +
    labs(linetype="p = q = N/2", colour="Factor") + xlab("Iteration") + scale_x_continuous(breaks = seq(0, 10, 2)) + ylim(c(0,1)) +
    labs(y = TeX("$\\sin\\Theta$ distance"))
  ggsave(g, file = paste(fig_path, sprintf("stat_err_converg_SNR%d_2.png", SNR/0.375), sep = ""),
         width = 8, height = 10, units = "cm")
}

#statistical error vs. 1/SNR
g <-ggplot(stat_err_results[stat_err_results$iteration==3,], aes(x = 1/SNR, y = est_err_mean, color = factor)) + 
  geom_errorbar(aes(ymin = est_err_mean - est_err_sd, ymax = est_err_mean + est_err_sd), width=.03) +
  geom_line(aes(x = 1/SNR, y = est_err_mean, color = factor, linetype = as.factor(dim)))  + ylim(c(0, 0.5)) +
  labs(linetype="p = q = N/2", colour="Factor", x = "1/SNR", y = TeX("$\\sin\\Theta$ distance"))

ggsave(g, file = paste(fig_path, "stat_err_SNR_2.png", sep = ""), width = 8, height = 10, units = "cm")


#combine wanted figures
results1_file <- "JisstPCA/simulations/TheoryValidation/JisstPCA_res.csv"
results2_file <- "JisstPCA/simulations/TheoryValidation/G-JisstPCA_res.csv"
fig_path <- "JisstPCA/simulations/TheoryValidation/fig/"

## validation of JisstPCA
ii <- 0
g <- list()
for(results_file in c(results1_file, results2_file)){
  ii <- ii + 1
  results = read.csv(results_file, header = TRUE)
  results_iter_summary  <- results %>% group_by(factor, dim, SNR, iteration)
  stat_err_results <- results_iter_summary %>% summarise(est_err_mean = mean(stat_err), est_err_sd = sd(stat_err))
  #statistical error vs. iteration
  SNR = 1.5
  g[[ii]] <- ggplot(stat_err_results[stat_err_results$SNR == SNR,], aes(x = iteration, y = est_err_mean, color = factor)) + 
    geom_errorbar(aes(ymin = est_err_mean - est_err_sd / sqrt(20) * qnorm(0.975), ymax = est_err_mean + est_err_sd / sqrt(20) * qnorm(0.975)), width=.2) +
    geom_line(aes(x = iteration, y = est_err_mean, color = factor, linetype = as.factor(dim))) +
    labs(linetype="p = q = N/2", colour="Factor") + xlab("Iteration") + scale_x_continuous(breaks = seq(0, 10, 2)) + ylim(c(0,0.6)) +
    labs(y = TeX("$\\sin\\Theta$ distance")) + theme(legend.position = "none") + theme(plot.margin=unit(c(0.7, rep(0.1,3)), "cm"))
  if(ii > 1){
    g[[ii]] <- g[[ii]] + ylab("") + theme(plot.margin=unit(c(0.7, 0.1, 0.1,-0.2), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  #statistical error vs. 1/SNR
  ii <- ii + 1
  g[[ii]] <-ggplot(stat_err_results[stat_err_results$iteration==3,], aes(x = 1/SNR, y = est_err_mean, color = factor)) + 
    geom_errorbar(aes(ymin = est_err_mean - est_err_sd, ymax = est_err_mean + est_err_sd), width=.03) +
    geom_line(aes(x = 1/SNR, y = est_err_mean, color = factor, linetype = as.factor(dim)))  + ylim(c(0, 0.6)) +
    labs(linetype="p = q = N/2", colour="Factor", x = "1/SNR", y = TeX("$\\sin\\Theta$ distance")) + ylab("") + 
    theme(plot.margin=unit(c(0.7, 0.1, 0.1,-0.2), "cm"), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if(ii < 4){
    g[[ii]] <- g[[ii]] + theme(legend.position = "none") 
  }
}

png(file = paste(fig_path, "Validation.png", sep = ""), width = 800, height = 300, units = "px")
print(grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]], widths = c(1.5,1.2,1.2,2), ncol = 4))
grid.text("JisstPCA", x = unit(0.26, "npc"), 
          y = unit(0.97, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
grid.text("G-JisstPCA", x = unit(0.65, "npc"), 
          y = unit(0.97, "npc"),gp = gpar(fontsize=15, fontfamily="Times New Roman"))
dev.off()




