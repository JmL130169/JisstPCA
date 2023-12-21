#summary for detailed results for generalized model
library(ggplot2)
library(cowplot)
library(patchwork)
library(pracma)
library(reshape)
library(latex2exp)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
fig_path <- "JisstPCA/simulations/ComparativeStudies/fig/generalModel_"

#read all results of generalized model (diagonal D)
#unstructured, diagonal D
singular_gap_setting <- 1
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/GeneralizedModel/diagonal_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results <- results[results$scenario == 4, ]
results <- results[results$seed <= 10, ]
results <- results[!results$orthogonal, ]
results$structure <- "unstructured"
results$d_setting <- "diagonal1"
results_combined <- results

#unstructured, diagonal D with more different eigenvalues
singular_gap_setting <- 2
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/GeneralizedModel/diagonal_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results <- results[results$seed <= 10, ]
results <- results[!results$orthogonal, ]
results$structure <- "unstructured"
results$d_setting <- "diagonal2"
results_combined <- rbind(results_combined, results)

#structured, diagonal D has small difference in eigenvalues
setting <- 1
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/StructuredSim/res%d.csv", setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results <- results[results$diagonal == 1, ]
results$structure <- "structured"
results$d_setting <- "diagonal1"
results_combined <- results_combined[ ,3:10]
results_combined <- rbind(results_combined, results[,2:9])

#structured, diagonal D, with larger difference in eigenvalues
setting <- 2
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/StructuredSim/res%d.csv", setting)
results = read.csv(res_filename, header = TRUE, sep = ",")
results <- results[!duplicated(results),]
results$structure <- "structured"
results$d_setting <- "diagonal2"
results_combined <- rbind(results_combined, results[,2:9])
results_combined <- results_combined[results_combined$SNR < 25, ]

results <- results_combined
results$method_type <- "ours_BIC_top"
results$method_type[results$method=="iHOSVD"|results$method=="iHOOI"] <- "baseline"
results$method_type[results$method=="JisstPCA0 (oracle)"|results$method=="G-JisstPCA0 (oracle)"|
                      results$method=="JisstPCA1 (oracle)"|results$method=="G-JisstPCA1 (oracle)"] <- "ours_oracle"
results$method_type[results$method=="JisstPCA0 (BIC)"|results$method=="G-JisstPCA0 (BIC)"|
                      results$method=="JisstPCA1 (BIC)"|results$method=="G-JisstPCA1 (BIC)"] <- "ours_BIC"
results$method[results$method == "JisstPCA0 (oracle)"] <- "JisstPCA (subtraction, oracle)"
results$method[results$method == "JisstPCA1 (oracle)"] <- "JisstPCA (projection, oracle)"
results$method[results$method == "G-JisstPCA0 (oracle)"] <- "G-JisstPCA (subtraction, oracle)"
results$method[results$method == "G-JisstPCA1 (oracle)"] <- "G-JisstPCA (projection, oracle)"

results$method[results$method == "JisstPCA0 (BIC)"] <- "JisstPCA (subtraction, BIC)"
results$method[results$method == "JisstPCA1 (BIC)"] <- "JisstPCA (projection, BIC)"
results$method[results$method == "G-JisstPCA0 (BIC)"] <- "G-JisstPCA (subtraction, BIC)"
results$method[results$method == "G-JisstPCA1 (BIC)"] <- "G-JisstPCA (projection, BIC)"

results$method[results$method == "JisstPCA0 (BIC top)"] <- "JisstPCA (subtraction, BIC top)"
results$method[results$method == "JisstPCA1 (BIC top)"] <- "JisstPCA (projection, BIC top)"
results$method[results$method == "G-JisstPCA0 (BIC top)"] <- "G-JisstPCA (subtraction, BIC top)"
results$method[results$method == "G-JisstPCA1 (BIC top)"] <- "G-JisstPCA (projection, BIC top)"

results$method[results$method == "iHOOI"] <- "iHOOI (oracle)"
results$method[results$method == "iHOSVD"] <- "iHOSVD (oracle)"

results_summary  <- results %>% group_by(SNR, structure, d_setting,method, factor.name, factor.number, method_type)
results_summary <- results_summary %>% summarise(err_mean = mean(err), err_sd = sd(err), nrep = n())
results_summary <- results_summary %>% unite("setting", 
                                             c(structure,d_setting), 
                                             remove = FALSE)

results_summary$method <- factor(results_summary$method, level = c("iHOSVD (oracle)", "iHOOI (oracle)", 
                                                                   "JisstPCA (subtraction, oracle)", 
                                                                   "G-JisstPCA (subtraction, oracle)", 
                                                                   "JisstPCA (projection, oracle)", 
                                                                   "G-JisstPCA (projection, oracle)",
                                                                   "JisstPCA (subtraction, BIC)", 
                                                                   "G-JisstPCA (subtraction, BIC)", 
                                                                   "JisstPCA (projection, BIC)", 
                                                                   "G-JisstPCA (projection, BIC)",
                                                                   "JisstPCA (subtraction, BIC top)", 
                                                                   "G-JisstPCA (subtraction, BIC top)", 
                                                                   "JisstPCA (projection, BIC top)", 
                                                                   "G-JisstPCA (projection, BIC top)"))

g <- list()
factor_number.labs <- c("factor 1", "factor 2"); names(factor_number.labs) <- c("1", "2")
setting_list <- c("unstructured_diagonal1", "unstructured_diagonal2", "structured_diagonal1", "structured_diagonal2" )
setting_names <- c("unstructured factors, generalized model setting 1", "unstructured factors, generalized model setting 2", 
                   "structured factors, generalized model setting 1", "structured factors, generalized model setting 2" )
for(method_setting in c("ours_oracle", "ours_BIC", "ours_BIC_top")){
    for(i in 1:4){
      g[[i]] <- ggplot(results_summary[(results_summary$method_type == method_setting|results_summary$method_type == "baseline")&
                                         results_summary$setting == setting_list[i],]) + 
        geom_line(aes(x = SNR, y = err_mean, color = method)) + 
        geom_errorbar(aes(x = SNR, ymin = err_mean - err_sd / sqrt(10) * qnorm(0.975), ymax = err_mean + err_sd / sqrt(nrep) * qnorm(0.975), color = method), width=0.6) +
        geom_point(aes(x = SNR, y = err_mean, color = method, shape = method), size = 1) + 
        facet_grid(factor.number~factor.name, switch = "y", scales = "free_y", labeller = labeller(factor.number = factor_number.labs)) + 
        theme(legend.position = "none",plot.title = element_text(size = 10, hjust = 0.5),
              plot.margin = margin(t = 0,  # Top margin
                                   r = 0.1,  # Right margin
                                   b = 0,  # Bottom margin
                                   l = 0,  # Left margin
                                   unit = "cm")) + 
        ylab(TeX("$\\sin\\Theta$ distance")) + xlab("SNR") + labs(title = setting_names[i])
    }
 
    g_combined <- ggarrange(g[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
             g[[2]] + theme(strip.text.y = element_blank(), axis.title.y = element_blank(), 
                         axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
              g[[3]] + theme(strip.text.x = element_blank()), 
             g[[4]] + theme(strip.text.x = element_blank(), strip.text.y = element_blank(), axis.title.y = element_blank()), 
              ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
    fig_filename <- paste(paste(fig_path, method_setting, sep = ""), ".png", sep = "")
    ggsave(g_combined, file = fig_filename, width = 20, height = 15, units = "cm", bg="white")
}

