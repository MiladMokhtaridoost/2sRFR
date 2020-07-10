library(Rmosek)
library(Matrix)

source("step1_L1_Linf.R")
source("step2_L1_decomposition.R")

DATA_PATH <- "./data"
RESULT_PATH <- "./results"
cohorts <- c("TCGA-ACC",  "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
             "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",  "TCGA-LIHC",
             "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",   "TCGA-PAAD",
             "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM",
             "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC",
             "TCGA-UCS",  "TCGA-UVM")
args <- commandArgs(trailingOnly = TRUE)
cohorts <- cohorts[as.numeric(args[[1]]) %% length(cohorts) + 1]

  load(sprintf("%s/TCGA-%s.RData", DATA_PATH, toupper(cohorts)))
  
  set.seed(6676)
####################preparation############  
  common_patident<-intersect(
    rownames(TCGA$mirna),rownames(TCGA$mrna))
  
  X<-TCGA$mirna[common_patident,]
  Y<-TCGA$mrna[common_patident,]
  
  X<-X[,which(colMeans(X>0)>=0.5)]
  X<-(scale(log2(X+1)))
  
  Y<-Y[,which(colMeans(Y>0)>=0.5)]
  Y<-(scale(log2(Y+1)))
  
#####################################
  
    n<-nrow(X) # Number of data instances (n)
    t<-ncol(Y) # Number of response variables (mRNAs)
    d <- ncol(X) # Number of predictor variables (miRNAs)  
    
    dimension <- list()
    dimension$n <- n
    dimension$t <- t
    dimension$d <- d  
  ################################
    
    regulatory_matrix <- step1_L1_Linf(dimension, X, Y)
    W <- regulatory_matrix$W
#####################################

  dimension$R_u <- 20 #upper bound for rank
  parameters <- list()
  parameters$maxiter <- 500 # maximum iteration in step 2 (decomposition)
  parameters$tol <- 1e-05 # convergence rate
  parameters$lam1 <- 3 # lambda_1
  parameters$lam2 <- 1 # lambda_2
  parameters$lam3 <- 4 # lambda_3
  parameters$lam4 <- 1.5 # lambda_4
  ###################decomposition######################
  
  decomp_out <- step2_L1_decomposition(dimension, parameters, W)
  Wx <- decomp_out$Wx
  Wy <- decomp_out$Wy
  R_star <- decomp_out$R_star
  iteration <- decomp_out$iteration
###############normalization###################
  colnames(Wy)<-colnames(Y,do.NULL = TRUE)
  rownames(Wx)<-colnames(X,do.NULL = TRUE)
  program_norms_normalized <- sqrt(colSums(Wx^2) * rowSums(Wy^2))
  sorted_indices_normalized <- sort(program_norms_normalized, decreasing = TRUE, index.return = TRUE)$ix
  
  summary <- list()
  summary$weights_normalized <- program_norms_normalized[sorted_indices_normalized]
  summary$Wx_normalized <- Wx[ ,sorted_indices_normalized]
  for (r in 1:ncol(Wx)) {
    summary$Wx_normalized[,r] <- summary$Wx_normalized[,r] / sqrt(sum(summary$Wx_normalized[,r]^2))
  }
  rownames(summary$Wx_normalized) <- rownames(Wx)
  
  summary$Wy_normalized <- Wy[sorted_indices_normalized,]
  for (r in 1:nrow(Wy)) {
    summary$Wy_normalized[r,] <- summary$Wy_normalized[r,] / sqrt(sum(summary$Wy_normalized[r,]^2))
  }
  colnames(summary$Wy_normalized) <- colnames(Wy)

#######save##########  
  summary$R_star <- R_star
  summary$W <- W
  summary$Wx <- Wx
  summary$Wy <- Wy
  summary$iteration <- iteration
  NRMSE <- mean(sqrt(colMeans((Wx %*% Wy - W)^2)) / apply(W, 2, sd))
  summary$obj <- decomp_out$obj_vals
  summary$NRMSE <- NRMSE
    
  save(summary, file = sprintf("%s/two-step_%s_summary.RData", RESULT_PATH, cohorts))
