MH_bc(Y, Y_zero, Cell_N, Gene_N, bc, Lambda, Pi, Theta, Nu, Sigma, step, burn, K){
  # Cell_N <- ncol(Y)
  # Gene_N <- nrow(Y)
  bc_sample <- matrix(0, ncol=K, nrow=Cell_N)
  ita <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Zexpbc <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bc <- sweep(Lambda, 2, exp(bc), FUN='*')
  Pi_bc <- sweep(Pi, 2, bc*Theta, FUN='+')
  Phi <- pnorm(Pi_bc)
  integrate_z <- apply(Y_zero*(1-Phi)+Phi*Lambda_bc^Y*exp(-Lambda_bc), 2, prod)
  poster_density <- exp(-(bc-Nu)^2/(2*Sigma^2))*integrate_z
  
  for(i in 1:burn){
    bc_propose <- bc + rnorm(Cell_N, mean=0, sd=step)
    Lambda_bc <- sweep(Lambda, 2, exp(bc_propose), FUN='*')
    Pi_bc <- sweep(Pi, 2, bc_propose*Theta, FUN='+')
    Phi <- pnorm(Pi_bc)
    integrate_z <- apply(Y_zero*(1-Phi)+Phi*Lambda^Y*exp(-Lambda), 2, prod)
    poster_density_propose <- exp(-(bc_propose-Nu)^2/(2*Sigma^2))*integrate_z
    ratio <- poster_density_propose/poster_density
    hidden <- runif(Cell_N)
    bc[ratio>hidden] <- bc_propose[ratio>hidden]
    poster_density[ratio>hidden] <- poster_density_propose[ratio>hidden]
  }
  for(i in 1:K){
    bc_propose <- bc + rnorm(Cell_N, mean=0, sd=step)
    Lambda_bc <- sweep(Lambda, 2, exp(bc_propose), FUN='*')
    Pi_bc <- sweep(Pi, 2, bc_propose*Theta, FUN='+')
    Phi <- pnorm(Pi_bc)
    integrate_z <- apply(Y_zero*(1-Phi)+Phi*Lambda^Y*exp(-Lambda), 2, prod)
    poster_density_propose <- exp(-(bc_propose-Nu)^2/(2*Sigma^2))*integrate_z
    ratio <- poster_density_propose/poster_density
    hidden <- runif(Cell_N)
    bc[ratio>hidden] <- bc_propose[ratio>hidden]
    poster_density[ratio>hidden] <- poster_density_propose[ratio>hidden]
    bc_sample[,i] <- bc
    Pi_bc <- sweep(Pi, 2, bc*Theta, FUN='+')
    Phi <- pnorm(Pi_bc)
    ita <- ita + Pi_bc + exp(-Pi_bc^2/2)/(sqrt(2*pi)*Phi)
    expbc <- exp(bc)
    Lambda_bc <- sweep(Lambda, 2, expbc, FUN='*')
    explambda <- exp(-Lambda_bc)
    z_temp <- Phi*explambda/(1+Phi*(explambda-1))
    Z <- Z + z_temp
    Zexpbc <- Zexpbc + sweep(z_temp, 2, expbc, FUN='*')
  }
  return(list(bc_sample=bc_sample, ita=ita, Z=Z, Zexpbc=Zexpbc))
}

fitDABEA <- function(Y, tot_read, bat_ind, bio_ind, gene_len, step=1, burn=1000, K=1000, max_iter=200, stop=0.001){
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
  Mu <- t(apply(Y, 1, FUN=function(x)aggregate(x, by=list(bio_ind), mean)[,2]))
  N_group <- length(unique(bio_ind))
  Mu1 <- matrix(0, nrow=Gene_N, ncol=N_group)
  bc <- rnrom(Cell_N, Nu, Sigma)

  Group_Matrix <- matrix(0, ncol=Cell_N, nrow=N_group)
  for(i in 1:N_group){
    Group_Matrix[i,]<-bio_ind==unique(bio_ind)[i]
  }
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  Pi <- sweep(sweep(matrix(Gamma, nrow=Gene_N, ncol=Cell_N), 2, Alpha*log(tot_read), FUN='+'), 1, Beta*log(gene_len), FUN='+')
  j <- 0
  while(j<=max_iter & max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), max(abs(Nu1-Nu)), max(abs(Sigma1-Sigma)), max(abs(Mu1-Mu)))>stop){
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
    
    Mu <- Y*Z%*%t(Group_Matrix)/sweep((sweep(Zexpbc, 2, tot_read, FUN='*')%*%t(Group_Matrix)), 1, gene_len, FUN='*')
    
    bc <- rowMeans(bc_sample)
    ita <- ita/K
    response <- as.vector(ita)
    bcvector <- rep(bc, each=Gene_N)
    LM_Results <- lm(response~logrc+loglg+bcvector)
    Gamma <- LM_Results$coefficients[[1]]
    Alpha <- LM_Results$coefficients[[2]]
    Beta <- LM_Results$coefficients[[3]]
    Theta <- LM_Results$coefficients[[4]]
    
    SumAll <- sum(bc_sample)
    SumEach <- c()
    SumSquare <- c()
    Batch_N <- c()
    Batch_Nu <- c()
    Batch_Nu0 <- rep(0, length(unique(bat_ind)))
    Batch_Sigma <- c()
    Batch_Sigma0 <- rep(1, length(unique(bat_ind)))
    for(batch in 1:length(unique(bat_ind))){
      batch_id <- unique(bat_ind)[batch]
      batch_bc <- bc_sample[bat_ind==batch_id,]
      SumEach <- c(SumEach, sum(batch_bc))
      SumSquare <- c(SumSquare, sum(batch_bc^2))
      Batch_N <- c(Batch_N, length(bat_ind==batch_id))
      Batch_Nu <- c(Batch_Nu, mean(batch_bc))
      Batch_Sigma <- c(Batch_Sigma, sd(batch_bc))
    }
    jj <- 0
    while(jj<=10 & max(max(abs(Batch_Sigma0-Batch_Sigma)), max(max(abs(Batch_Nu0-Batch_Nu))))>=stop){
      Batch_Sigma0 <- Batch_Sigma
      Batch_Nu0 <- Batch_Nu
      temp <- sum(Batch_Sigma*Batch_N)
      Batch_Nu <- (SumEach-Batch_Sigma*Batch_N*SumAll/temp)/Batch_N
      Batch_Sigma <- (SumSquare-2*Batch_Nu*SumEach+Batch_N*Batch_Nu^2)/Batch_N
      jj <- jj + 1
    }
    Nu <- rep(0, Cell_N)
    Sigma <- rep(0, Cell_N)
    for(batch in 1:length(unique(bat_ind))){
      batch_id <- unique(bat_ind)[batch]
      Nu[bat_ind==batch_id] <- Batch_Nu[batch]
      Sigma[bat_ind==batch_id] <- Batch_Sigma[batch]
    }
    j <- j + 1
  }
  return(list(Mu=Mu, Alpha=Alpha, Beta=Beta, Gamma=Gamma, Theta=Theta, Nu=Nu, Sigma=Sigma, bc=bc))
}