library(ggplot2)
library(cowplot)
library(patchwork)
library(pracma)
library(reshape)
library(latex2exp)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggpubr)

U_filename <- "JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorU_q80_N20_seed2.csv"
VW_filename <- "JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorVW_q80_N20_seed2.csv"
VW_weighted_filename <- "JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorVW_weighted_q80_N20_seed2.csv"
u_results <- read.csv(U_filename, header = TRUE, sep = ",")
VW_results <- read.csv(VW_filename, header = TRUE, sep = ",")
VW_weighted_results <- read.csv(VW_filename, header = TRUE, sep = ",")
#u_results <- read.csv(U_filename, header = FALSE, sep = ",")
#colnames(u_results) <- c("method", "FactorNumber", "RowInd", "Val")
#VW_results <- read.csv(VW_filename, header = FALSE, sep = ",")
#colnames(VW_results) <- c("method", "FactorType", "FactorName", "FactorNumber", 
#                          "RowInd", "ColInd", "Val", "OrderedRowInd", "OrderedColInd")

#VW_results_combined <- rbind(VW_results[VW_results$method == "iHOSVD"|
#                                          VW_results$method == "iHOOI",], 
#                             VW_weighted_results[VW_results$method != "iHOSVD"&
#                                                   VW_results$method != "iHOOI",])


VW_results <- VW_results %>% unite("FactorType", c(FactorName,FactorNumber), 
                                   remove = FALSE)
VW_results$method[VW_results$method == "truth"] <- "Truth"
VW_results$method[VW_results$method == "iHOSVD"] <- "iHOSVD (oracle)"
VW_results$method[VW_results$method == "iHOOI"] <- "iHOOI (oracle)"
VW_results$method <- factor(VW_results$method, levels = c("Truth", "JisstPCA (oracle)", "JisstPCA (BIC)", 
                                                          "G-JisstPCA (oracle)", "G-JisstPCA (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)"))
method_list <- c("Truth", "JisstPCA (BIC)", 
                 "G-JisstPCA (BIC)", "iHOSVD (oracle)", "iHOOI (oracle)")
VW_results <- VW_results[VW_results$method%in%method_list, ]
colors <- colorRampPalette(c("lightyellow", "dodgerblue1", "dodgerblue2", "dodgerblue3"))

plot_v1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="V",], 
                 aes(x = RowInd, y = ColInd, fill = Val)) + 
  geom_tile() + scale_fill_gradient(low = "lightyellow", high = "dodgerblue3") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + theme(legend.position = "none") +
  #theme(legend.key.height= unit(0.5, 'cm'),
  #      legend.key.width= unit(0.5, 'cm')) + 
  #theme(legend.title = element_text(size=8)) + 
  #theme(legend.text  = element_text(size=5)) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(method~., switch = "y") + 
  ggtitle(TeX("$V_{1}V_{1}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(0,-0.4,0,-0.4), "cm"))
plot_v1
plot_v2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="V",], 
                 aes(x = RowInd, y = ColInd, fill = Val)) + 
  geom_tile() + scale_fill_gradient(low = "lightyellow", high = "dodgerblue3") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.position = "none") +
  theme(strip.text.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(method~., switch = "y") + 
  ggtitle(TeX("$V_{2}V_{2}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(0,-0.4,0,-0.4), "cm"))
plot_v2

plot_w1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="W",], 
                 aes(x = RowInd, y = ColInd, fill = Val)) + 
  geom_tile() + scale_fill_gradient(low = "lightyellow", high = "dodgerblue3") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.position = "none") +
  theme(strip.text.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(method~., switch = "y") + 
  ggtitle(TeX("$W_{1}W_{1}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(0,-0.4,0,-0.4), "cm"))
plot_w1

plot_w2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="W",], 
                 aes(x = RowInd, y = ColInd, fill = Val)) + 
  geom_tile() + scale_fill_gradient(low = "lightyellow", high = "dodgerblue3") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.position = "none") +
  theme(strip.text.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(method~., switch = "y") + 
  ggtitle(TeX("$W_{2}W_{2}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.margin=unit(c(0,-0.4,0,-0.4), "cm"))
plot_w2


png(file = "JisstPCA/simulations/ComparativeStudies/NetworkSim/networkVisResFull.png", width = 400, height = 500, units = "px")
print(grid.arrange(plot_v1, plot_v2, plot_w1, plot_w2,
                   widths = c(1.4, 1.2, 1.2, 1.2), nrow = 1, ncol = 4))
dev.off()

