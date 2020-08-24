library(Rmosek)
step2_L1_decomposition <- function(dimension, parameters, W) {
 
  R <- dimension$R_u 
  n <- dimension$n 
  t <- dimension$t 
  d <- dimension$d  
  
  maxiter <- parameters$maxiter
  tol <- parameters$tol 
  lam1 <- parameters$lam1
  lam2 <- parameters$lam2
  lam3 <- parameters$lam3
  lam4 <- parameters$lam4
################################  
 obj_vals <- matrix(0,nrow = maxiter, ncol=1)
 Wx <- Matrix(runif((R * d),min = -1,max = 1), ncol = R, nrow = d)

 ##################################  
 Wy1 <- Matrix(0, ncol=t, nrow=R)
 Wy2 <- Matrix(0, ncol=t, nrow=R)
 Wx1 <- Matrix(0, nrow=d, ncol=R)
 Wx2 <- Matrix(0, nrow=d, ncol=R)

 rowg <- matrix(c(1,-1),nrow = 1)
 A0<- bdiag(replicate((d*t), rowg, simplify = FALSE))

 obj_vals[1] <- 99999999
#####################start iterative algorithm#######################
iter=1
 while (iter <= maxiter) {
  iter = iter+1
  ################Wy LP#########################
    print(iter)
    A1 <- bdiag(replicate(t, Wx, simplify = FALSE))
    A1 <- cbind(A1, -A1)

    AWy <- cbind(A1,A0)
    CWy <- cbind(matrix(lam3, ncol = (2*t*R)), matrix(1, ncol = (2*d*t)))
    C0Wy <- lam1*(sum(abs(Wx))) + lam2 * sqrt(sum(Wx^2))
    lbWy <- matrix(0,ncol=(2*t*(d+R)))
    ubWy <- matrix(Inf, ncol = (2*t*(d+R))) 

    ################################Solve for Wy#############################
    problem_Wy <- list()
    problem_Wy$sense <- "min"
    problem_Wy$c <- CWy
    problem_Wy$c0 <- C0Wy
    problem_Wy$A <- AWy
    problem_Wy$bc <- rbind(blc = c(W) , buc = c(W))
    problem_Wy$bx <- rbind(blx = lbWy, bux = ubWy) 
    problem_Wy$qobj <- list(i = c(1:R),
                         j = c(1:R),
                         v = rep(lam4, R))
    
    opts <- list()  
    problem_Wy$iparam = list(LOG=0)
    result <- mosek(problem_Wy, opts)
    Wy1 <- matrix(result$sol$itr$xx[1:(R*t)],ncol=t)
    Wy2 <- matrix(result$sol$itr$xx[((R*t)+1):(2*R*t)],ncol=t)
    Wy <- Wy1-Wy2
    
    ################################Wx LP##############################
    A2 <- bdiag(replicate(d, t(Wy), simplify = FALSE))
    A2 <- cbind(A2,-A2)
    
    AWx <- cbind(A2,A0)
    CWx <- cbind(matrix(lam1, ncol = (2*d*R)), matrix(1, ncol = (2*d*t)))
    C0Wx <- lam3*(sum(abs(Wy)) + lam4 * sqrt(sum(Wy^2)))
    lbWx <- matrix(0,ncol=(2*d*(t+R)))
    ubWx <- matrix(Inf, ncol = (2*d*(t+R))) 
    
################################Solve for Wx#############################
    problem_Wx <- list()
    problem_Wx$sense <- "min"
    problem_Wx$c <- t(CWx)
    problem_Wx$c0 <- C0Wx
    problem_Wx$A <- AWx
    problem_Wx$bc <- rbind(blc = c(t(W)) , buc = c(t(W)))
    problem_Wx$bx <- rbind(blx = lbWx, bux = ubWx) 
    problem_Wx$qobj <- list(i = c(1:d),
                         j = c(1:d),
                         v = rep(lam2, d))
    
    opts <- list()
    problem_Wx$iparam = list(LOG=0)
    result <- mosek(problem_Wx, opts)
    Wx1 <- matrix(result$sol$itr$xx[1:(R*d)],nrow=R)
    Wx2 <- matrix(result$sol$itr$xx[((R*d)+1):(2*R*d)],nrow=R)
    Wx <- t(Wx1)-t(Wx2)
    ##########check stoping creteria & rank condition##########

    obj_vals[iter] <- sum(abs(W -  Wx %*% Wy))+ lam1*sum(abs(Wx)) + lam2*sqrt(sum(Wx^2)) +lam3*sum(abs(Wy)) + lam4*sqrt(sum(Wy^2))
 
    dif <- (abs(obj_vals[iter]- obj_vals[iter-1]))/obj_vals[iter-1] 
    if (dif < tol || is.nan(dif)) {
      if (((qr(Wy)$rank != R) || (qr(Wx)$rank != R))) {
        R <- R-1
        
        Wx <- Matrix(runif((R * d),min = -1,max = 1), ncol = R, nrow = d)
        Wy1 <- Matrix(0, ncol=t, nrow=R)
        Wy2 <- Matrix(0, ncol=t, nrow=R)
        Wx1 <- Matrix(0, nrow=d, ncol=R)
        Wx2 <- Matrix(0, nrow=d, ncol=R)
      }
      else  {
        
        break
      }
    }
}
 decomp_out <- list(Wx=Wx, Wy=Wy, R_star=R, obj_vals=obj_vals, iteration=iter)
}
 
