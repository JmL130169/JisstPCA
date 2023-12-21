library(ggplot2)
library(cowplot)
library(patchwork)
library(pracma)
library(reshape)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(tidyr)
fig_path <- "/Applications/Files/research/JisstPCA/JisstPCA_new/fig/"

#first collect results from all experiments
#unstructured, scalar d, 1/2 singular gap
singular_gap_setting <- 1
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/UnstructuredSim/scalar_d_res1.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE)
results <- results[!duplicated(results),]
results <- results[results$scenario == 4, ]
results$structure <- "unstructured"
results$d_setting <- "scalar"
#orthogonal with siginificant singular gap and non-orthogonal
results$orthogonal_setting <- ifelse(results$orthogonal, "orthogonal2", "non-orthogonal")
results_combined <- results

#unstructured, scalar d, small singular gap
singular_gap_setting <- 2
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/UnstructuredSim/scalar_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE)
results <- results[!duplicated(results),]
results <- results[results$scenario == 4, ]
results$structure <- "unstructured"
results$d_setting <- "scalar"
#orthogonal with little or no singular gap
results$orthogonal_setting <- "orthogonal1"
results_combined <- rbind(results_combined, results)

#unstructured, diagonal D
singular_gap_setting <- 1
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/GeneralizedModel/diagonal_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results <- results[results$scenario == 4, ]
results <- results[results$seed <= 10, ]
results$structure <- "unstructured"
results$d_setting <- "diagonal1"
results$orthogonal_setting <- ifelse(results$orthogonal, "orthogonal1", "non-orthogonal")
results_combined <- rbind(results_combined, results)

#unstructured, diagonal D with more different eigenvalues
singular_gap_setting <- 2
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/GeneralizedModel/diagonal_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results <- results[results$seed <= 10, ]
results$structure <- "unstructured"
results$d_setting <- "diagonal2"
results$orthogonal_setting <- ifelse(results$orthogonal, "orthogonal1", "non-orthogonal")
results_combined <- rbind(results_combined, results)

#structured, diagonal D has small difference in eigenvalues
res_filename <- "JisstPCA/simulations/ComparativeStudies/StructuredSim/res.csv"
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results$structure <- "structured"
results$d_setting <- ifelse(results$diagonal, "diagonal1", "scalar")
results$orthogonal_setting <- "non-orthogonal"
results_combined <- results_combined[ ,3:11]
results_combined <- rbind(results_combined, results[,2:10])
results_combined$method[results_combined$method == "iHOOI"] <- "iHOOI (oracle)"
results_combined$method[results_combined$method == "iHOSVD"] <- "iHOSVD (oracle)"

results_summary  <- results_combined %>% group_by(SNR, structure, d_setting, orthogonal_setting, method, factor.name, factor.number)
results_summary <- results_summary %>% summarise(err_mean = mean(err), err_sd = sd(err), nrep = n())
results_summary <- results_summary %>% unite("setting", 
                                             c(structure,orthogonal_setting,d_setting), 
                                             remove = FALSE)

#plotted setting
setting_list <- c("unstructured_non-orthogonal_scalar", "unstructured_orthogonal1_scalar", 
                  "structured_non-orthogonal_scalar", "unstructured_non-orthogonal_diagonal1")
results_summary <- results_summary[results_summary$setting %in% setting_list, ]


subset <- which(results_summary$method=="iHOOI (oracle)"|results_summary$method=="iHOSVD (oracle)"|
                  results_summary$method == "G-JisstPCA0 (BIC)"|results_summary$method == "JisstPCA0 (BIC)")
results_sub <- results_summary[subset,]
results_sub$method[results_sub$method=="JisstPCA0 (BIC)"] <- "JisstPCA (BIC)"; 
results_sub$method[results_sub$method=="G-JisstPCA0 (BIC)"] <- "G-JisstPCA (BIC)"
results_sub$method <- factor(results_sub$method, level = c("iHOSVD (oracle)", "iHOOI (oracle)", "JisstPCA (BIC)", "G-JisstPCA (BIC)"))
factor_number.labs <- c("factor 1", "factor 2"); names(factor_number.labs) <- c("1", "2")

results_sub <- results_sub[results_sub$SNR < 25, ]
setting_names <- c("unstructured and non-orthogonal factors", "unstructured and orthogonal factors", 
                   "structured and non-orthogonal factors", "unstructured and non-orthogonal factors, generalized model")
g <- list()
for(i in 1 : 4){
  g[[i]] <- ggplot(results_sub[results_sub$setting == setting_list[i],]) + 
    geom_line(aes(x = SNR, y = err_mean, color = method)) + 
    geom_errorbar(aes(x = SNR, ymin = err_mean - err_sd / sqrt(10) * qnorm(0.975), ymax = err_mean + err_sd / sqrt(nrep) * qnorm(0.975), color = method), 
                  width=0.6) +
    geom_point(aes(x = SNR, y = err_mean, shape = method, color = method), size = 1) + 
    facet_grid(factor.number~factor.name, switch = "y", scales = "free_y", labeller = labeller(factor.number = factor_number.labs)) + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 10, hjust = 0.5),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0.1,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0,  # Left margin
                               unit = "cm")) + 
    ylab(TeX("$\\sin\\Theta$ distance")) + xlab("SNR") + labs(title = setting_names[i], )
}
g_combined <- ggarrange(g[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
                        g[[2]] + theme(strip.text.y = element_blank(), axis.title.y = element_blank(), 
                                       axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
                        g[[3]] + theme(strip.text.x = element_blank()), 
                        g[[4]] + theme(strip.text.x = element_blank(), strip.text.y = element_blank(), axis.title.y = element_blank()), 
                        ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
fig_filename <- paste(fig_path, "CombinedComp_main.png", sep = "")
ggsave(g_combined, file = fig_filename, width = 20, height = 15, units = "cm", bg="white")

