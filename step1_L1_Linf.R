step1_L1_Linf <- function(dimension, X, Y) {

W1 <- Matrix(0, ncol=1, nrow=dimension$d)
W2 <- Matrix(0, ncol=1, nrow=dimension$d)
W <- Matrix(0, ncol=dimension$t, nrow=dimension$d)
AW<- Matrix(0, ncol=(2*dimension$d), nrow=(2*dimension$n))

###########L1/inf regression###########
AW[1:dimension$n , 1:d] <- X
AW[1:dimension$n , (d+1):(2*d)] <- -X
AW[(dimension$n+1):(2*dimension$n) , 1:dimension$d] <- -X
AW[(dimension$n+1):(2*dimension$n) , (dimension$d+1):(2*dimension$d)] <- X
AW <- cbind(AW,1)

CW <- cbind(matrix(0, ncol = (2*dimension$d)),1)  

lbW <- matrix(0,ncol=(2*dimension$d+1))
ubW <- matrix(Inf, ncol = (2*dimension$d+1))  

################################Solve for W#############################
for (i in 1:dimension$t){  
  problem_W <- list()
  problem_W$sense <- "min"
  problem_W$c <- t(CW)
  problem_W$A <- AW
  problem_W$bc <- rbind(blc = cbind(t(Y[ ,i]),-t(Y[ ,i])) , buc = rep(Inf,(2*dimension$n)))
  problem_W$bx <- rbind(blx = lbW, bux = ubW) 
  opts <- list() 
  problem_W$iparam = list(LOG=0)
  result <- mosek(problem_W, opts)
  W1 <- result$sol$itr$xx[1:dimension$d]
  W2 <- result$sol$itr$xx[(dimension$d+1):(2*dimension$d)]
  W[,i] <- W1-W2
} 
W <- as.matrix(W)
regulatory_matrix <- list(W = W)
}
