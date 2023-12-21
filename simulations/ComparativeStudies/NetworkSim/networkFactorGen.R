library(ggplot2)
library(cowplot)
library(patchwork)
library(pracma)
library(reshape)

r = c(3, 2); p = 80
within_prob = matrix(c(0.8, 0.6, 0.7, 0.5), 2, 2, byrow = TRUE)
off_prob <- 0.3
B_x2 = matrix(rep(off_prob, r[2]^2), r[2], r[2]) + diag(rep(within_prob[1,2] - off_prob, r[2]))
B_x1 = matrix(rep(off_prob, r[1]^2), r[1], r[1]) + diag(rep(within_prob[1,1] - off_prob, r[1]))
B_y2 = matrix(rep(off_prob, r[2]^2), r[2], r[2]) + diag(rep(within_prob[2,2] - off_prob, r[2]))
B_y1 = matrix(rep(off_prob, r[1]^2), r[1], r[1]) + diag(rep(within_prob[2,1] - off_prob, r[1]))
for(q in c(50, 80)){
  set.seed(10)
  theta_x1 = matrix(0, nrow = p, ncol = r[1])
  theta_x2 = matrix(0, nrow = p, ncol = r[2])
  theta_y1 = matrix(0, nrow = q, ncol = r[1])
  theta_y2 = matrix(0, nrow = q, ncol = r[2])
  clu_x1 = c(rep(1, 0.4*p), rep(2, 0.3*p), rep(3, 0.3*p))
  #clu_x1 = clu_x1[sample(1:p, replace=FALSE)]
  clu_x2 = c(rep(1, 0.3*p), rep(2, 0.2*p), rep(1, 0.2*p), rep(2, 0.3*p))
  #clu_x2 = clu_x2[sample(1:p, replace=FALSE)]
  clu_y1 = c(rep(1, 0.4*q), rep(2, 0.4*q), rep(3, 0.2*q))
  #clu_y1 = clu_y1[sample(1:q, replace=FALSE)]
  clu_y2 = c(rep(1, 0.3*q), rep(2, 0.4*q), rep(1, 0.3*q))
  #clu_y2 = clu_y2[sample(1:q, replace=FALSE)]
  
  for (i in 1:p){
    theta_x1[i, clu_x1[i]] = 1
    theta_x2[i, clu_x2[i]] = 1
  }
  for (i in 1:q){
    theta_y1[i, clu_y1[i]] = 1
    theta_y2[i, clu_y2[i]] = 1
  }
  
  # first layer of X
  px1 = theta_x1%*%B_x1%*%t(theta_x1)
  # second layer of X
  px2 = theta_x2%*%B_x2%*%t(theta_x2)
  # first layer of Y
  py1 = theta_y1%*%B_y1%*%t(theta_y1) 
  # second layer of Y
  py2 = theta_y2%*%B_y2%*%t(theta_y2)
  write.csv(px1, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_px1_q%d.csv", q), row.names = F)
  write.csv(px2, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_px2_q%d.csv", q), row.names = F)
  write.csv(py1, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_py1_q%d.csv", q), row.names = F)
  write.csv(py2, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_py2_q%d.csv", q), row.names = F)
  write.csv(clu_x1, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_blocknumber_V1_q%d.csv", q), row.names = F)
  write.csv(clu_x2, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_blocknumber_V2_q%d.csv", q), row.names = F)
  write.csv(clu_y1, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_blocknumber_W1_q%d.csv", q), row.names = F)
  write.csv(clu_y2, file = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_blocknumber_W2_q%d.csv", q), row.names = F)
}


