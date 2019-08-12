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
  # integrate_z <- apply(Y_zero*(1-Phi)+Phi*Lambda_bc^Y*exp(-Lambda_bc), 2, prod)
  # poster_density <- exp(-(bc-Nu)^2/(2*Sigma^2))*integrate_z
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
    #integrate_z <- apply((Y_zero*(1-Phi_propose)+Phi_propose*(Lambda_bc_propose/Lambda_bc)^Y*exp(-Lambda_bc_propose))/(Y_zero*(1-Phi)+Phi*exp(-Lambda_bc)), 2, prod)
    #poster_density_propose <- exp(-(bc_propose-Nu)^2/(2*Sigma^2))*integrate_z
    #ratio <- poster_density_propose/poster_density
    #ratio <- exp((-(bc_propose-Nu)^2+(bc-Nu)^2)/(2*Sigma^2))*integrate_z
    log_ratio <- sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0))
    ratio <- exp(log_ratio)
    hidden <- runif(Cell_N)
    accept <- ratio>hidden
    bc[accept] <- bc_propose[accept]
    # Lambda_bc[,accept] <- Lambda_bc_propose[,accept]
    # Pi_bc[,accept] <- Pi_bc_propose[,accept]
    # Phi[,accept] <- Phi_propose[,accept]
    log_dense[accept] <- log_dense_propose[accept]
    #poster_density[ratio>hidden] <- poster_density_propose[ratio>hidden]
  }
  for(i in 1:K){
    bc_propose <- bc + rnorm(Cell_N, mean=0, sd=step)
    Lambda_bc_propose <- sweep(Lambda, 2, exp(bc_propose), FUN='*')
    Pi_bc_propose <- sweep(Pi, 2, bc_propose*Theta, FUN='+')
    Phi_propose <- pnorm(Pi_bc_propose)
    log_int_z[Y_zero] <- log(1+Phi_propose[Y_zero]*(exp(-Lambda_bc_propose[Y_zero])-1))
    log_int_z[!Y_zero] <- log(Phi_propose[!Y_zero])+Y[!Y_zero]*log(Lambda_bc_propose[!Y_zero])-Lambda_bc_propose[!Y_zero]
    log_dense_propose <- -(bc_propose-Nu)^2/(2*Sigma^2)+colSums(log_int_z)
    #integrate_z <- apply((Y_zero*(1-Phi_propose)+Phi_propose*(Lambda_bc_propose/Lambda_bc)^Y*exp(-Lambda_bc_propose))/(Y_zero*(1-Phi)+Phi*exp(-Lambda_bc)), 2, prod)
    #poster_density_propose <- exp(-(bc_propose-Nu)^2/(2*Sigma^2))*integrate_z
    #ratio <- poster_density_propose/poster_density
    log_ratio <- sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0))
    ratio <- exp(log_ratio)
    hidden <- runif(Cell_N)
    accept <- ratio>hidden
    bc[accept] <- bc_propose[accept]
    Lambda_bc[,accept] <- Lambda_bc_propose[,accept]
    Pi_bc[,accept] <- Pi_bc_propose[,accept]
    Phi[,accept] <- Phi_propose[,accept]
    #poster_density[ratio>hidden] <- poster_density_propose[ratio>hidden]
    log_dense[accept] <- log_dense_propose[accept]
    bc_sample[,i] <- bc
    #Pi_bc <- sweep(Pi, 2, bc*Theta, FUN='+')
    #Phi <- pnorm(Pi_bc)
    ita <- ita + Pi_bc + exp(-Pi_bc^2/2)/(sqrt(2*pi)*Phi)
    expbc <- exp(bc)
    #Lambda_bc <- sweep(Lambda, 2, expbc, FUN='*')
    explambda <- exp(-Lambda_bc)
    z_temp <- Phi*explambda/(1+Phi*(explambda-1))
    Z <- Z + z_temp
    Zexpbc <- Zexpbc + sweep(z_temp, 2, expbc, FUN='*')
  }
  return(list(bc_sample=bc_sample, ita=ita, Z=Z, Zexpbc=Zexpbc))
}

fitDABEA <- function(Y, tot_read, bat_ind, bio_ind, gene_len, step=0.01, burn=50, burn_start=500, K=500, max_iter=200, stop=0.005){
 gene_len <- gene_len/mean(gene_len)
 tot_read <- tot_read/mean(tot_read)
  Y_zero <- Y <= 0
  Cell_N <- ncol(Y)
  Gene_N <- nrow(Y)
  logrc <- rep(log(tot_read), each=Gene_N)
  loglg <- rep(log(gene_len), Cell_N)
  
  #start point
  Gamma <- 0
  Gamma1 <- 1
  Alpha <- 0
  Alpha1 <- 1
  Beta <- 0
  Beta1 <- 1
  Theta <- 0
  Theta1 <- 1
  Nu <- rep(0, Cell_N)
  Nu1 <- rep(1, Cell_N)
  Sigma <- rep(1, Cell_N)
  Sigma1 <- rep(1, Cell_N)
  N_group <- length(unique(bio_ind))
  Mu <- matrix(apply(Y, 1, FUN=function(x)aggregate(x, by=list(bio_ind), mean)[,2]), ncol = N_group) + 0.01
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
  
  j <- 1  
  
  while(j<=max_iter & max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma), abs(Mu1-Mu))>stop){
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
    # SumAll <- sum(bc_sample)
    # SumEach <- c()
    # SumSquare <- c()
    # Batch_N <- c()
    # Batch_Nu <- c()
    # Batch_Nu0 <- rep(0, length(unique(bat_ind)))
    # Batch_Sigma <- c()
    # Batch_Sigma0 <- rep(1, length(unique(bat_ind)))
    # for(batch in 1:length(unique(bat_ind))){
    #   batch_id <- unique(bat_ind)[batch]
    #   batch_bc <- bc_sample[bat_ind==batch_id,]
    #   SumEach <- c(SumEach, sum(batch_bc))
    #   SumSquare <- c(SumSquare, sum(batch_bc^2))
    #   Batch_N <- c(Batch_N, length(bat_ind==batch_id))
    #   Batch_Nu <- c(Batch_Nu, mean(batch_bc))
    #   Batch_Sigma <- c(Batch_Sigma, sd(batch_bc))
    # }
    # jj <- 0
    # while(jj<=100 & max(abs(Batch_Sigma0-Batch_Sigma), abs(Batch_Nu0-Batch_Nu))>=stop){
    #   Batch_Sigma0 <- Batch_Sigma
    #   Batch_Nu0 <- Batch_Nu
    #   temp <- sum(Batch_Sigma*Batch_N)
    #   Batch_Nu <- (SumEach-Batch_Sigma*Batch_N*SumAll/temp)/Batch_N
    #   Batch_Sigma <- (SumSquare-2*Batch_Nu*SumEach+Batch_N*Batch_Nu^2)/Batch_N
    #   jj <- jj + 1
    # }
    # # Nu <- rep(0, Cell_N)
    # # Sigma <- rep(0, Cell_N)
    # for(batch in 1:length(unique(bat_ind))){
    #   batch_id <- unique(bat_ind)[batch]
    #   Nu[bat_ind==batch_id] <- Batch_Nu[batch]
    #   Sigma[bat_ind==batch_id] <- Batch_Sigma[batch]
    # }
    
    j <- j + 1
  }
  
  meanNu <- mean(Nu)
  bc = bc - meanNu
  Nu = Nu - meanNu
  Mu = Mu * exp(meanNu)
  Gamma = Gamma + Theta*meanNu
  
  return(list(Mu=Mu, Alpha=Alpha, Beta=Beta, Gamma=Gamma, Theta=Theta, Nu=Nu, Sigma=Sigma, bc=bc))
}