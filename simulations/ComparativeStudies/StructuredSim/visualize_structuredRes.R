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

U_filename <- "JisstPCA/simulations/ComparativeStudies/StructuredSims/factorU.csv"
VW_filename <- "JisstPCA/simulations/ComparativeStudies/StructuredSims/factorVW.csv"
u_results <- read.csv(U_filename, header = TRUE, sep = ",")
VW_results <- read.csv(VW_filename, header = TRUE, sep = ",")

VW_results <- VW_results %>% unite("FactorType", c(FactorName,FactorNumber), 
                                                   remove = FALSE)
VW_results$method <- factor(VW_results$method, levels = c("truth", "JisstPCA", "G-JisstPCA", "iHOSVD", "iHOOI"))

colors <- colorRampPalette(c("lightyellow", "dodgerblue1", "dodgerblue2", "dodgerblue3", "darkblue"))

plot_v1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="V",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_v2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="V",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_w1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="W",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_w2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="W",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

u_results$method <- factor(u_results$method, levels = c("truth", "JisstPCA", "G-JisstPCA", "iHOSVD", "iHOOI"))
u_res1 <- u_results[u_results$FactorNumber==1, c(1, 3, 4)]; u_res2 <- u_results[u_results$FactorNumber==2, c(1, 3, 4)];
colnames(u_res1)[3] <- "u1_val"; colnames(u_res2)[3] <- "u2_val";
u_res <- merge(u_res1, u_res2, by = c("method", "RowInd"), all = TRUE)
u_cluster <- read.csv("/Applications/Files/research/JisstPCA/JisstPCA_new/structure/u_cluster.csv", header = TRUE, sep = ",")
u_res$cluster <- u_cluster[u_res$RowInd, 1]
plot_u = ggplot(u_res, aes(x=u1_val, y=u2_val, color=cluster)) + geom_point(shape=1) + ggtitle("Scatter plot of u") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
  facet_grid(method~., switch = "y") + xlab(TeX("$u_1$")) + ylab(TeX("$u_2$")) + 
  theme(aspect.ratio = 1.2, plot.margin=unit(c(0,-0.6,0,-0.85), "cm")) 
plot_u
png(file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/VisResFull.png", width = 550, height = 500, units = "px")
print(grid.arrange(plot_v1, plot_v2, plot_w1, plot_w2, plot_u,
                   widths = c(0.7, 0.6, 0.6, 0.6, 0.85), nrow = 1, ncol = 5))
dev.off()


##plot results from main methods
VW_results <- VW_results[VW_results$method == "truth"|VW_results$method == "JisstPCA"|VW_results$method == "iHOOI", ]
u_res <- u_res[u_res$method == "truth"|u_res$method == "JisstPCA"|u_res$method == "iHOOI", ]
plot_v1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="V",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_v2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="V",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_w1 = ggplot(VW_results[VW_results$FactorNumber==1&VW_results$FactorName=="W",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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

plot_w2 = ggplot(VW_results[VW_results$FactorNumber==2&VW_results$FactorName=="W",], aes(x = RowInd, y = ColInd, fill = abs(Val))) + 
  geom_tile() + scale_fill_gradientn(colors = colors(100)) + theme_bw() +
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
u
plot_u = ggplot(u_res, aes(x=u1_val, y=u2_val, color=cluster)) + geom_point(shape=1) + ggtitle("Scatter plot of u") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") +
  facet_grid(method~., switch = "y") + xlab(TeX("$u_1$")) + ylab(TeX("$u_2$")) + 
  theme(aspect.ratio = 1.2, plot.margin=unit(c(0,-0.6,0,-0.85), "cm")) 
plot_u
png(file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/VisResMain.png", width = 830, height = 500, units = "px")
print(grid.arrange(plot_v1, plot_v2, plot_w1, plot_w2, plot_u,
                   widths = c(0.71, 0.64, 0.64, 0.64, 0.75), nrow = 1, ncol = 5))
dev.off()



