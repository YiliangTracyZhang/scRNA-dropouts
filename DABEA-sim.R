
#####################
### fit function ###
####################

## Metropolis hasting algorithm

MH_bc <- function(Y, Y_zero, Cell_N, Gene_N, bc, Lambda, Pi, Theta, Nu, Sigma, step, burn, K){
  # Cell_N <- ncol(Y)
  # Gene_N <- nrow(Y)
  bc_sample <- matrix(0, ncol=K, nrow=Cell_N)
  ita <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Zexpbc <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bc <- sweep(Lambda, 2, exp(bc), FUN='*')
  Pi_bc <- sweep(Pi, 2, bc*Theta, FUN='+')
  Phi <- pnorm(Pi_bc)
  log_int_z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bc[Y_zero])-1))
  log_int_z[!Y_zero] <- log(Phi[!Y_zero])+Y[!Y_zero]*log(Lambda_bc[!Y_zero])-Lambda_bc[!Y_zero]
  log_dense <- -(bc-Nu)^2/(2*Sigma^2)+colSums(log_int_z)
  for(i in 1:burn){
    bc_propose <- bc + rnorm(Cell_N, mean=0, sd=step)
    Lambda_bc <- sweep(Lambda, 2, exp(bc_propose), FUN='*')
    Pi_bc <- sweep(Pi, 2, bc_propose*Theta, FUN='+')
    Phi <- pnorm(Pi_bc)
    log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bc[Y_zero])-1))
    log_int_z[!Y_zero] <- log(Phi[!Y_zero])+Y[!Y_zero]*log(Lambda_bc[!Y_zero])-Lambda_bc[!Y_zero]
    log_dense_propose <- -(bc_propose-Nu)^2/(2*Sigma^2)+colSums(log_int_z)
    log_ratio <- sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0))
    ratio <- exp(log_ratio)
    hidden <- runif(Cell_N)
    accept <- ratio>hidden
    bc[accept] <- bc_propose[accept]
    log_dense[accept] <- log_dense_propose[accept]
  }
  for(i in 1:K){
    bc_propose <- bc + rnorm(Cell_N, mean=0, sd=step)
    Lambda_bc_propose <- sweep(Lambda, 2, exp(bc_propose), FUN='*')
    Pi_bc_propose <- sweep(Pi, 2, bc_propose*Theta, FUN='+')
    Phi_propose <- pnorm(Pi_bc_propose)
    log_int_z[Y_zero] <- log(1+Phi_propose[Y_zero]*(exp(-Lambda_bc_propose[Y_zero])-1))
    log_int_z[!Y_zero] <- log(Phi_propose[!Y_zero])+Y[!Y_zero]*log(Lambda_bc_propose[!Y_zero])-Lambda_bc_propose[!Y_zero]
    log_dense_propose <- -(bc_propose-Nu)^2/(2*Sigma^2)+colSums(log_int_z)
    log_ratio <- sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0))
    ratio <- exp(log_ratio)
    hidden <- runif(Cell_N)
    accept <- ratio>hidden
    bc[accept] <- bc_propose[accept]
    Lambda_bc[,accept] <- Lambda_bc_propose[,accept]
    Pi_bc[,accept] <- Pi_bc_propose[,accept]
    Phi[,accept] <- Phi_propose[,accept]
    log_dense[accept] <- log_dense_propose[accept]
    bc_sample[,i] <- bc
    expbc <- exp(bc)
    explambda <- exp(-Lambda_bc)
    z_temp <- Phi*explambda/(1+Phi*(explambda-1))
    z_temp[!Y_zero] <- 1
    Z <- Z + z_temp
    ita <- ita + Pi_bc + z_temp*exp(-Pi_bc^2/2)/(sqrt(2*pi)*Phi) - (1-z_temp)*exp(-Pi_bc^2/2)/(sqrt(2*pi)*(1-Phi))  
    Zexpbc <- Zexpbc + sweep(z_temp, 2, expbc, FUN='*')
  }
  return(list(bc_sample=bc_sample, ita=ita, Z=Z, Zexpbc=Zexpbc))
}

## main function 

fitDABEA <- function(Y, tot_read, bat_ind, bio_ind, gene_len, step=0.1, burn=50, burn_start=500, K=500, max_iter=200, stop=0.01){
 # gene_len <- gene_len/mean(gene_len)
 # tot_read <- tot_read/mean(tot_read)
  Y_zero <- Y <= 0
  Cell_N <- ncol(Y)
  Gene_N <- nrow(Y)
  logrc <- rep(log(tot_read), each=Gene_N)
  loglg <- rep(log(gene_len), Cell_N)
  
  #start point
  Gamma <- -0.5
  Gamma1 <- 1
  Alpha <- 0.1
  Alpha1 <- 0
  Beta <- 0.1
  Beta1 <- 0
  Theta <- 1
  Theta1 <- 0
  Nu <- rep(0, Cell_N)
  Nu1 <- rep(1, Cell_N)
  Sigma <- rep(1, Cell_N)
  Sigma1 <- rep(1, Cell_N)
  N_group <- length(unique(bio_ind))
  # for real data
  # Mu <- matrix(apply(Y, 1, FUN=function(x)aggregate(x, by=list(bio_ind), mean)[,2]), ncol = N_group) + 0.01 ##
  # for simulated data
  Mu <- matrix(0,nrow=Gene_N,ncol=N_group)
  for(i in 1:N_group){
          Mu[,i]<-matrix(Y[,bio_ind==unique(bio_ind)[i]],nrow = Gene_N) %*% (1/tot_read[bio_ind==unique(bio_ind)[i]] )* (1/gene_len)
          Mu[,i] <- Mu[,i] / (length(bio_ind[bio_ind == unique(bio_ind)[i]]))
  }
  Mu1 <- matrix(0, nrow=Gene_N, ncol=N_group)
  bc <- rnorm(Cell_N, Nu, Sigma)

  Group_Matrix <- matrix(0, ncol=Cell_N, nrow=N_group)
  for(i in 1:N_group){
    Group_Matrix[i,]<-bio_ind==unique(bio_ind)[i]
  }
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  Pi <- sweep(sweep(matrix(Gamma, nrow=Gene_N, ncol=Cell_N), 2, Alpha*log(tot_read), FUN='+'), 1, Beta*log(gene_len), FUN='+')
  
  bat_indexes <- unique(bat_ind)
  N_bat <- length(bat_indexes)
  
  # E step
  MH <- MH_bc(Y, Y_zero, Cell_N, Gene_N, bc, Lambda, Pi, Theta, Nu, Sigma, step, burn_start, K)
  bc_sample <- MH$bc_sample
  ita <- MH$ita
  ita <- ita/K
  Z <- MH$Z
  Z <- Z/K
  Zexpbc <- MH$Zexpbc
  Zexpbc <- Zexpbc/K
  
       #####################
  # observe the traceplot
  trace <- MH$bc_sample
  trace_1 <- trace[1,]
  p <- plot(x=1:length(trace_1), y=trace_1, type = 'n', main = paste0('trace plot for 1 iteration'), ylab = 'value', xlab = 'sampling times') + lines(x = 1:length(trace_1), y = trace_1)  
  print(p)
        ######################
  
  #M step
  Gamma1 <- Gamma
  Alpha1 <- Alpha
  Beta1 <- Beta
  Theta1 <- Theta
  Nu1 <- Nu
  Sigma1 <- Sigma
  Mu1 <- Mu
  
  Mu <- (Y*Z)%*%t(Group_Matrix)/sweep((sweep(Zexpbc, 2, tot_read, FUN='*')%*%t(Group_Matrix)), 1, gene_len, FUN='*')
  
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*') #### 08/16 newly added
  Pi <- sweep(sweep(matrix(Gamma, nrow=Gene_N, ncol=Cell_N), 2, Alpha*log(tot_read), FUN='+'), 1, Beta*log(gene_len), FUN='+') ### 08/16 newly added
  
  bc <- rowMeans(bc_sample)
  # ita <- ita/K
  response <- as.vector(ita)
  bcvector <- rep(bc, each=Gene_N)
  LM_Results <- lm(response~logrc+loglg+bcvector)
  Gamma <- LM_Results$coefficients[[1]]
  Alpha <- LM_Results$coefficients[[2]]
  Beta <- LM_Results$coefficients[[3]]
  Theta <- LM_Results$coefficients[[4]]
  
  for(batch in 1:N_bat){
          batch_id <- bat_indexes[batch]
          batch_bc <- bc_sample[bat_ind==batch_id,]
          Nu[bat_ind==batch_id] <- mean(batch_bc)
          Sigma[bat_ind==batch_id] <- sd(batch_bc)
  }
  # meanNu <- mean(Nu)
  # bc = bc - meanNu
  # Nu = Nu - meanNu
  # Mu = Mu * exp(meanNu)
  # Gamma = Gamma + Theta*meanNu
  j <- 1  
  
  print(c(date(),paste0('the', j , 'iteration')))
  print(c('error', max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma))))
  print(c('Mu error', max(abs(Mu1-Mu))))
  print(Mu[1:10,])
  print(Alpha)
  print(Beta)
  print(Gamma)
  print(Theta)
  print(unique(Nu))
  print(unique(Sigma))
  
  while(j<=max_iter & max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma))>stop){
    # E step
    MH <- MH_bc(Y, Y_zero, Cell_N, Gene_N, bc, Lambda, Pi, Theta, Nu, Sigma, step, burn, K)
    bc_sample <- MH$bc_sample
    ita <- MH$ita
    Z <- MH$Z
    Zexpbc <- MH$Zexpbc
    #M step
    Gamma1 <- Gamma
    Alpha1 <- Alpha
    Beta1 <- Beta
    Theta1 <- Theta
    Nu1 <- Nu
    Sigma1 <- Sigma
    Mu1 <- Mu
    
    Mu <- (Y*Z)%*%t(Group_Matrix)/sweep((sweep(Zexpbc, 2, tot_read, FUN='*')%*%t(Group_Matrix)), 1, gene_len, FUN='*')
    
    Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*') ####
    Pi <- sweep(sweep(matrix(Gamma, nrow=Gene_N, ncol=Cell_N), 2, Alpha*log(tot_read), FUN='+'), 1, Beta*log(gene_len), FUN='+') ###
    
    bc <- rowMeans(bc_sample)
    ita <- ita/K
    response <- as.vector(ita)
    bcvector <- rep(bc, each=Gene_N)
    LM_Results <- lm(response~logrc+loglg+bcvector)
    Gamma <- LM_Results$coefficients[[1]]
    Alpha <- LM_Results$coefficients[[2]]
    Beta <- LM_Results$coefficients[[3]]
    Theta <- LM_Results$coefficients[[4]]

    for(batch in 1:N_bat){
      batch_id <- bat_indexes[batch]
      batch_bc <- bc_sample[bat_ind==batch_id,]
      Nu[bat_ind==batch_id] <- mean(batch_bc)
      Sigma[bat_ind==batch_id] <- sd(batch_bc)
    }
    # meanNu <- mean(Nu)
    # bc = bc - meanNu
    # Nu = Nu - meanNu
    # Mu = Mu * exp(meanNu)
    # Gamma = Gamma + Theta*meanNu
    
    print(c(date(),paste0('the', j , 'iteration')))
    print(c('error', max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma))))
    print(c('Mu error', max(abs(Mu1-Mu))))
    print(Mu[1:10,])
    print(Alpha)
    print(Beta)
    print(Gamma)
    print(Theta)
    print(unique(Nu))
    print(unique(Sigma))
    
    j <- j + 1
  }
  
  print(paste0(j,' iterations in total'))
  print(c('error', max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma))))
  print(c('Mu error', max(abs(Mu1-Mu))))
  print(Mu[1:10,])
  print(Alpha)
  print(Beta)
  print(Gamma)
  print(Theta)
  print(unique(Nu))
  print(unique(Sigma))
  
  #####dropout imputation
  Y_impute <- Y + Lambda[Y_zero]
  Z_matrix <- Z
  ##### batch effect sampling (stored for quality control)
  sampling <- bc_sample
  
  return(list(Mu=Mu, 
              Alpha=Alpha, Beta=Beta, Gamma=Gamma, Theta=Theta, Nu=Nu, Sigma=Sigma, 
              Y_impute = Y_impute, Z_matrix = Z_matrix,
              bc=bc, sampling = sampling))
}