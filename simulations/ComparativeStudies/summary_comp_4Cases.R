#summary for detailed results for 4 dimensionality
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
fig_path <- "JisstPCA/simulations/ComparativeStudies/fig/"

#unstructured, scalar d, 1/2 singular gap; 
#show dimension vs sample size effects and orthogonal with singular gap comparison results
singular_gap_setting <- 1
res_filename <- sprintf("JisstPCA/simulations/ComparativeStudies/UnstructuredSim/scalar_d_res%d.csv", singular_gap_setting)
results = read.csv(res_filename, header = TRUE)
results <- results[!duplicated(results),]
results <- results[results$SNR < 25,]
results$structure <- "unstructured"
results$d_setting <- "scalar"
#orthogonal with siginificant singular gap and non-orthogonal
results$orthogonal_setting <- ifelse(results$orthogonal, "orthogonal2", "non-orthogonal")
results$method_type <- ""
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

results$method[results$method == "iHOOI"] <- "iHOOI (oracle)"
results$method[results$method == "iHOSVD"] <- "iHOSVD (oracle)"

results_summary  <- results %>% group_by(scenario, SNR, structure, d_setting, orthogonal_setting, method, 
                                         method_type, factor.name, factor.number)
results_summary <- results_summary %>% summarise(err_mean = mean(err), err_sd = sd(err), nrep = n())

results_summary$method <- factor(results_summary$method, level = c("iHOSVD (oracle)", "iHOOI (oracle)", 
                                                           "JisstPCA (subtraction, oracle)", 
                                                           "G-JisstPCA (subtraction, oracle)", 
                                                           "JisstPCA (projection, oracle)", 
                                                           "G-JisstPCA (projection, oracle)",
                                                           "JisstPCA (subtraction, BIC)", 
                                                           "G-JisstPCA (subtraction, BIC)", 
                                                           "JisstPCA (projection, BIC)", 
                                                           "G-JisstPCA (projection, BIC)"))
scenario.labs <- c("case 1", "case 2", "case 3", "case 4");names(scenario.labs) <- c("1", "2", "3", "4")
g <- list()
for(orth_setting in c("non-orthogonal", "orthogonal2")){
  for(method_setting in c("ours_oracle", "ours_BIC")){
    for(i in 1:2){
      g[[i]] <- ggplot(results_summary[results_summary$orthogonal_setting == orth_setting&
                                         (results_summary$method_type == method_setting|results_summary$method_type == "baseline")&
                                        results_summary$factor.number==i,]) + 
        geom_line(aes(x = SNR, y = err_mean, color = method)) + 
        geom_errorbar(aes(x = SNR, ymin = err_mean - err_sd / sqrt(10) * qnorm(0.975), ymax = err_mean + err_sd / sqrt(nrep) * qnorm(0.975), color = method), width=0.6) +
        geom_point(aes(x = SNR, y = err_mean, color = method, shape = method), size = 1) + 
        facet_grid(scenario~factor.name, switch = "y", labeller = labeller(scenario = scenario.labs)) + 
        theme(legend.position = "none",plot.title = element_text(size = 10, hjust = 0.5),
         plot.margin = margin(t = 0,  # Top margin
                               r = 0.1,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0,  # Left margin
                               unit = "cm")) + 
       ylab(TeX("$\\sin\\Theta$ distance")) + xlab("SNR") + labs(title = sprintf("Factor %d", i))
    }
    g_combined <- ggarrange(g[[1]],
                            g[[2]] + theme(strip.text.y = element_blank(), axis.title.y = element_blank()),
                            ncol=2, common.legend = TRUE, legend="bottom")
    
    #g_combined <- ggarrange(g[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
    #         g[[2]] + theme(strip.text.y = element_blank(), axis.title.y = element_blank(), 
    #                     axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
    #          g[[3]] + theme(strip.text.x = element_blank()), 
    #         g[[4]] + theme(strip.text.x = element_blank(), strip.text.y = element_blank(), axis.title.y = element_blank()), 
    #          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
    fig_filename <- paste(paste(paste(fig_path, orth_setting, sep = ""), method_setting, sep = "_"), ".png", sep = "")
    ggsave(g_combined, file = fig_filename, width = 20, height = 15, units = "cm", bg="white")
  }
}

