n = 200
p = c(0.3, 0.4, 0.3)
mu1 = c(0, 5, 15)
mu2 = c(0, 15, 5)
sigma = rep(1, 3)

set.seed(123)
ind_1 = runif(n, 0, 1)
u_1 = rep(0, n)
u_2 = rep(0, n)
cluster = c()
for (i in 1:n){
  u_1[i] = rnorm(1, mu1[1], sigma[1])*(ind_1[i] < p[1]) + rnorm(1, mu1[2], sigma[2])*(ind_1[i] > p[1] && ind_1[i] < (p[1] + p[2])) + rnorm(1, mu1[3], sigma[3])*(ind_1[i] > (p[1] + p[2]))
  u_2[i] = rnorm(1, mu2[1], sigma[1])*(ind_1[i] < p[1]) + rnorm(1, mu2[2], sigma[2])*(ind_1[i] > p[1] && ind_1[i] < (p[1] + p[2])) + rnorm(1, mu2[3], sigma[3])*(ind_1[i] > (p[1] + p[2]))
  if (ind_1[i] < p[1]){
    cluster[i] = "group 1"
  } else if (ind_1[i] > p[1] && ind_1[i] < (p[1] + p[2])){
    cluster[i] = "group 2"
  } else{
    cluster[i] = "group 3"
  }
}

u_1 = u_1/norm(as.matrix(u_1))
u_2 = u_2/norm(as.matrix(u_2))

gs_u = gramSchmidt(cbind(u_1, u_2), tol = 0.001)
u = gs_u$Q
u_1 = u[, 1]
u_2 = u[, 2]
u = data.frame(u)

u = data.frame(cbind(u_1, u_2))
plot_u_t = ggplot(u, aes(x=u_1, y=u_2, color=cluster)) + geom_point(shape=1) + ggtitle("Scatter plot of true u") + theme(plot.title = element_text(hjust = 0.5)) + theme(aspect.ratio=1) + theme(legend.position = "none")
plot_u_t

# factors of X
px = 50
r = c(3, 2)
px_1 = c(20, 20, 10)
px_2 = c(25, 25)
spx_1 = cumsum(px_1)
spx_2 = cumsum(px_2)
V1 = matrix(0, px, r[1])
V2 = matrix(0, px, r[2])

set.seed(2023)
V1[1:px_1[1], 1] = rnorm(px_1[1], 0, 1)
for (i in 1:(r[1]-1)){
  V1[(spx_1[i]+1):(spx_1[i+1]), i+1] = rnorm(px_1[i+1], 0, 1)
}
V1 = orth(V1)

# one version of V2: block gaussian
# V2[1:px_2[1], 1] = rnorm(px_2[1], 5, 25)
# for (i in 1:(r[2]-1)){
#   V2[(spx_2[i]+1):(spx_2[i+1]), i+1] = rnorm(px_2[i+1], 5, 25)
# }
# V2 = orth(V2)

# another version of V2: star graph
V_V2 = matrix(0, px, px)
L_V = matrix(0, 25, 25)
L_V[1, 1] = 24
for (i in 2:25){
  L_V[1, i] = -1
  L_V[i, 1] = -1
  L_V[i, i] = 1
}
V_V2[1:25, 1:25] = L_V
V_V2[26:50, 26:50] = L_V

v2_svd = svd(V_V2)
V2 = v2_svd$u[, 1:r[2]]

# factors of Y
py = 50
r = c(3, 2)
py_1 = c(20, 15, 15)
py_2 = c(30, 20)
spy_1 = cumsum(py_1)
spy_2 = cumsum(py_2)
W1 = matrix(0, py, r[1])
W2 = matrix(0, py, r[2])

set.seed(2023)
W1[1:py_1[1], 1] = rnorm(py_1[1], 0, 1)
for (i in 1:(r[1]-1)){
  W1[(spy_1[i]+1):(spy_1[i+1]), i+1] = rnorm(py_1[i+1], 0, 1)
}
W1 = orth(W1)

# one version of V2: block gaussian
# V2[1:px_2[1], 1] = rnorm(px_2[1], 5, 25)
# for (i in 1:(r[2]-1)){
#   V2[(spx_2[i]+1):(spx_2[i+1]), i+1] = rnorm(px_2[i+1], 5, 25)
# }
# V2 = orth(V2)

# another version of V2: star graph
W_W2 = matrix(0, py, py)
L_W1 = matrix(0, 30, 30)
L_W2 = matrix(0, 20, 20)
L_W1[1, 1] = 29
L_W2[1, 1] = 19
for (i in 2:30){
  L_W1[1, i] = -1
  L_W1[i, 1] = -1
  L_W1[i, i] = 1
}
for (i in 2:20){
  L_W2[1, i] = -1
  L_W2[i, 1] = -1
  L_W2[i, i] = 1
}
W_W2[1:30, 1:30] = L_W1
W_W2[31:50, 31:50] = L_W2

w2_svd = svd(W_W2)
W2 = w2_svd$u[, 1:r[2]]
write.csv(u_1, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/u1.csv", row.names = F)
write.csv(u_2, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/u2.csv", row.names = F)
write.csv(V1, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/V1.csv", row.names = F)
write.csv(V2, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/V2.csv", row.names = F)
write.csv(W1, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/W1.csv", row.names = F)
write.csv(W2, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/W2.csv", row.names = F)
write.csv(cluster, file = "JisstPCA/simulations/ComparativeStudies/StructuredSims/u_cluster.csv", row.names = F)
# plot and export of V1V1' and V2V2'
# V1V1'
V_V_1 = abs(V1%*%t(V1))
colnames(V_V_1) <- paste(1:50)
rownames(V_V_1) <- paste(1:50)

# Transform the matrix in long format
df_1 <- melt(V_V_1)
colnames(df_1) <- c("x", "y", "value")
plot_v1 = ggplot(df_1, aes(x = x, y = y, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + coord_fixed() +
  theme(axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 5), 
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) + 
  theme(legend.title = element_text(size=8)) + 
  theme(legend.text  = element_text(size=5)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   
  ggtitle(TeX("$V_{1}V_{1}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5))

# V2V2'
V_V_2 = abs(V2%*%t(V2))
colnames(V_V_2) <- paste(1:50)
rownames(V_V_2) <- paste(1:50)

# Transform the matrix in long format
df_2 <- melt(V_V_2)
colnames(df_2) <- c("x", "y", "value")
plot_v2 = ggplot(df_2, aes(x = x, y = y, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + coord_fixed() +
  theme(axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 5), 
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) + 
  theme(legend.title = element_text(size=8)) + 
  theme(legend.text  = element_text(size=5)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   
  ggtitle(TeX("$V_{2}V_{2}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5))

# W1W1'
W_W_1 = abs(W1%*%t(W1))
colnames(W_W_1) <- paste(1:50)
rownames(W_W_1) <- paste(1:50)

# Transform the matrix in long format
df_3 <- melt(W_W_1)
colnames(df_3) <- c("x", "y", "value")
plot_w1 = ggplot(df_3, aes(x = x, y = y, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + coord_fixed() +
  theme(axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 5), 
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) + 
  theme(legend.title = element_text(size=8)) + 
  theme(legend.text  = element_text(size=5)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   
  ggtitle(TeX("$W_{1}W_{1}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5))

# W2W2'
W_W_2 = abs(W2%*%t(W2))
colnames(W_W_2) <- paste(1:50)
rownames(W_W_2) <- paste(1:50)

# Transform the matrix in long format
df_4 <- melt(W_W_2)
colnames(df_4) <- c("x", "y", "value")
plot_w2 = ggplot(df_4, aes(x = x, y = y, fill = value)) + geom_tile() + scale_fill_gradient(low = "white", high = "red") + coord_fixed() +
  theme(axis.title.x = element_text(size = 8), 
        axis.text.x = element_text(size = 5), 
        axis.title.y = element_text(size = 8), 
        axis.text.y = element_text(size = 5), 
        title = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(aspect.ratio=1) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm')) + 
  theme(legend.title = element_text(size=8)) + 
  theme(legend.text  = element_text(size=5)) +
  theme(plot.margin = margin(0, 0, 0, 0)) +   
  ggtitle(TeX("$W_{2}W_{2}^{\\prime}$")) + theme(plot.title = element_text(hjust = 0.5))

plot_case1 = plot_v1 + plot_v2 + plot_w1 + plot_w2
plot_case1
