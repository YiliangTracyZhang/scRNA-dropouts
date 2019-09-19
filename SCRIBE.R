####################################
### Metropolis hasting algorithm ###
####################################

MH_bg <- function(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg, Lambda, Pi, Theta, Nu, Sigma, step, burn, K){
  bg_sum <- matrix(0, ncol=N_batch, nrow=Gene_N)
  bgsq_sum <- rep(0, N_batch)
  ita <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Zexpbg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Lambda_bg_propose <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Pi_bg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  Pi_bg_propose <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg[,batches]), FUN='*')
    Pi_bg[, bat_ind == batch_name[batches]] <- sweep(Pi[, bat_ind == batch_name[batches]], 1, bg[,batches]*Theta[batches], FUN='+')
  }
  Phi <- pnorm(Pi_bg)
  log_int_z <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bg[Y_zero])-1))
  log_int_z[!Y_zero] <- log(Phi[!Y_zero])+Y[!Y_zero]*log(Lambda_bg[!Y_zero])-Lambda_bg[!Y_zero]
  log_dense <- matrix(0, ncol=N_batch, nrow=Gene_N)
  log_dense_propose <- matrix(0, ncol=N_batch, nrow=Gene_N)
  for(batches in 1:N_batch){
    log_dense[, batches] <- -(bg[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
  }
  for(i in 1:burn){
    bg_propose <- bg + matrix(rnorm(N_batch*Gene_N, mean=0, sd=step), ncol=N_batch)
    for(batches in 1:N_batch){
      Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_propose[,batches]), FUN='*')
      Pi_bg[, bat_ind == batch_name[batches]] <- sweep(Pi[, bat_ind == batch_name[batches]], 1, bg_propose[,batches]*Theta[batches], FUN='+')
    }
    Phi <- pnorm(Pi_bg)
    log_int_z[Y_zero] <- log(1+Phi[Y_zero]*(exp(-Lambda_bg[Y_zero])-1))
    log_int_z[!Y_zero] <- log(Phi[!Y_zero])+Y[!Y_zero]*log(Lambda_bg[!Y_zero])-Lambda_bg[!Y_zero]
    for(batches in 1:N_batch){
      log_dense_propose[, batches] <- -(bg_propose[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
    }
    log_ratio <- matrix(sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0)), ncol=N_batch)
    ratio <- exp(log_ratio)
    hidden <- matrix(runif(N_batch*Gene_N), ncol=N_batch)
    accept <- ratio>hidden
    bg[accept] <- bg_propose[accept]
    log_dense[accept] <- log_dense_propose[accept]
  }
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg[,batches]), FUN='*')
    Pi_bg[, bat_ind == batch_name[batches]] <- sweep(Pi[, bat_ind == batch_name[batches]], 1, bg[,batches]*Theta[batches], FUN='+')
  }
  Phi <- pnorm(Pi_bg)
  for(i in 1:K){
    bg_propose <- bg + matrix(rnorm(N_batch*Gene_N, mean=0, sd=step), ncol=N_batch)
    for(batches in 1:N_batch){
      Lambda_bg_propose[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_propose[,batches]), FUN='*')
      Pi_bg_propose[, bat_ind == batch_name[batches]] <- sweep(Pi[, bat_ind == batch_name[batches]], 1, bg_propose[,batches]*Theta[batches], FUN='+')
    }
    Phi_propose <- pnorm(Pi_bg_propose)
    log_int_z[Y_zero] <- log(1+Phi_propose[Y_zero]*(exp(-Lambda_bg_propose[Y_zero])-1))
    log_int_z[!Y_zero] <- log(Phi_propose[!Y_zero])+Y[!Y_zero]*log(Lambda_bg_propose[!Y_zero])-Lambda_bg_propose[!Y_zero]
    for(batches in 1:N_batch){
      log_dense_propose[, batches] <- -(bg_propose[,batches]-Nu[batches])^2/(2*Sigma[batches]^2)+rowSums(log_int_z[,bat_ind == batch_name[batches]])
    }
    log_ratio <- matrix(sapply(log_dense_propose - log_dense, FUN = function(x)min(x,0)), ncol=N_batch)
    ratio <- exp(log_ratio)
    hidden <- matrix(runif(N_batch*Gene_N), ncol=N_batch)
    accept <- ratio>hidden
    bg[accept] <- bg_propose[accept]
    for(batches in 1:N_batch){
      Lambda_bg[accept[,batches], bat_ind == batch_name[batches]] <- Lambda_bg_propose[accept[,batches],bat_ind == batch_name[batches]]
      Pi_bg[accept[,batches], bat_ind == batch_name[batches]] <- Pi_bg_propose[accept[,batches], bat_ind == batch_name[batches]]
      Phi[accept[,batches], bat_ind == batch_name[batches]] <- Phi_propose[accept[,batches], bat_ind == batch_name[batches]]
    }
    log_dense[accept] <- log_dense_propose[accept]
    bg_sum <- bg_sum + bg
    expbg <- exp(bg)
    explambda <- exp(-Lambda_bg)
    z_temp <- Phi*explambda/(1+Phi*(explambda-1))
    z_temp[!Y_zero] <- 1
    Z <- Z + z_temp
    ita <- ita + Pi_bg + z_temp*exp(-Pi_bg^2/2)/(sqrt(2*pi)*Phi) - (1-z_temp)*exp(-Pi_bg^2/2)/(sqrt(2*pi)*(1-Phi))  
    for(batches in 1:N_batch){
      Zexpbg[,bat_ind == batch_name[batches]] <- Zexpbg[,bat_ind == batch_name[batches]] + sweep(z_temp[, bat_ind == batch_name[batches]], 1, expbg[,batches], FUN='*')
      bgsq_sum[batches] <- bgsq_sum[batches] + sum(bg[,batches]^2)
    }
  }
  return(list(bg_sum=bg_sum, bgsq_sum=bgsq_sum, ita=ita, Zexpbc=Zexpbc, Z=Z))
}

#####################
### main function ###
#####################

fitSCRIBE <- function(Y, bat_ind, bio_ind, step=0.1, burn=50, burn_start=500, K=500, max_iter=100, stop=0.001){
  # gene_len <- gene_len/mean(gene_len)
  # tot_read <- tot_read/mean(tot_read)
  Y_zero <- Y <= 0
  Cell_N <- ncol(Y)
  Gene_N <- nrow(Y)
  gene_len <- rowMeans(Y)/mean(Y)
  # logrc <- rep(log(tot_read), each=Gene_N)
  # loglg <- rep(log(gene_len), Cell_N)
  group_name <- unique(bio_ind)
  N_group <- length(group_name)
  batch_name = unique(bat_ind)
  N_batch <- length(batch_name)
  
  tot_read <- rep(0, Cell_N)
  for(batches in 1:N_batch){
    tot_read[bat_ind==batch_name[batches]] <- colMeans(Y[,bat_ind==batch_name[batches]])/mean(Y[,bat_ind==batch_name[batches]])
  }
  
  #start point
  Gamma <- rep(0, N_batch)
  Alpha <- rep(0, N_batch)
  Beta <- rep(0, N_batch)
  Theta <- rep(0, N_batch)
  Nu <- rep(0, N_batch)
  Sigma <- rep(0.5, N_batch)
  Mu <- matrix(0,nrow=Gene_N,ncol=N_group)
  for(groups in 1:N_group){
    Mu[,groups] <- rowMeans(sweep(Y[,bio_ind==group_name[groups]], 2, 1/tot_read[bio_ind==group_name[groups]], FUN='*'))/(gene_len)
  }
  bg <- matrix(0, ncol=N_batch, nrow=Gene_N)
  for(batches in 1:N_batch){
    bg[,batches] <- rnorm(Gene_N, mean=Nu[batches], sd=Sigma[batches])
  }
  
  Group_Matrix <- matrix(0, ncol=Cell_N, nrow=N_group)
  for(groups in 1:N_group){
    Group_Matrix[groups,]<-bio_ind==group_name[groups]
  }
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  Pi <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi, 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]], FUN="+")), 1, Beta[batches]*log(gene_len), FUN="+")
  }
  
  # E step
  MH <- MH_bg(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg, Lambda, Pi, Theta, Nu, Sigma, step, burn_start, K)
  bg_av <- MH$bg_sum/K
  bgsq_av <- MH$bgsq_sum/(K*Gene_N)
  ita <- MH$ita/K
  Zexpbc <- MH$Zexpbc
  
  #M step
  Gamma1 <- Gamma
  Alpha1 <- Alpha
  Beta1 <- Beta
  Theta1 <- Theta
  Nu1 <- Nu
  Sigma1 <- Sigma
  Mu1 <- Mu
  
  for(groups in 1:N_group){
    Mu[,groups] <- K*rowSums(Y[,bio_ind==group_name[groups]])/(gene_len*rowSums(sweep(Zexpbg[,bio_ind==group_name[groups]], 2, tot_read[bio_ind==group_name[groups]], FUN="*")))
  }
  
  for(batches in 1:N_batch){
    batch_ita <- ita[,bat_ind == batch_name[batches]]
    response <- as.vector(batch_ita)
    bgvector <- rep(bg[,batches], sum(bat_ind == batch_name[batches]))
    logrc <- rep(log(tot_read[bat_ind == batch_name[batches]]), each=Gene_N)
    loglg <- rep(log(gene_len), sum(bat_ind == batch_name[batches]))
    LM_Results <- lm(response~logrc+loglg+bgvector)
    Gamma[batches] <- LM_Results$coefficients[[1]]
    Alpha[batches] <- LM_Results$coefficients[[2]]
    Beta[batches] <- LM_Results$coefficients[[3]]
    Theta[batches] <- LM_Results$coefficients[[4]]
  }
  
  for(batches in 1:N_batch){
    Nu[batch] <- mean(bg_av[,batches])
    Sigma[batch] <- bgsq_av[batches] - Nu[batch]^2
  }
  Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
  for(batches in 1:N_batch){
    Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi, 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]], FUN="+")), 1, Beta[batches]*log(gene_len), FUN="+")
  }
  j <- 1
  while(j<=max_iter & max(abs(Alpha1-Alpha), abs(Beta1-Beta), abs(Theta1-Theta), abs(Nu1-Nu), abs(Sigma1-Sigma),abs(Mu-Mu1))>stop){
    # E step
    MH <- MH_bg(Y, Y_zero, Cell_N, Gene_N, N_batch, batch_name, bat_ind, bg_av, Lambda, Pi, Theta, Nu, Sigma, step, burn, K)
    bg_av <- MH$bg_sum/K
    bgsq_av <- MH$bgsq_sum/(K*Gene_N)
    ita <- MH$ita/K
    Zexpbc <- MH$Zexpbc
    
    #M step
    Gamma1 <- Gamma
    Alpha1 <- Alpha
    Beta1 <- Beta
    Theta1 <- Theta
    Nu1 <- Nu
    Sigma1 <- Sigma
    Mu1 <- Mu
    
    for(groups in 1:N_group){
      Mu[,groups] <- K*rowSums(Y[,bio_ind==group_name[groups]])/(gene_len*rowSums(sweep(Zexpbg[,bio_ind==group_name[groups]], 2, tot_read[bio_ind==group_name[groups]], FUN="*")))
    }
    
    for(batches in 1:N_batch){
      batch_ita <- ita[,bat_ind == batch_name[batches]]
      response <- as.vector(batch_ita)
      bgvector <- rep(bg[,batches], sum(bat_ind == batch_name[batches]))
      logrc <- rep(log(tot_read[bat_ind == batch_name[batches]]), each=Gene_N)
      loglg <- rep(log(gene_len), sum(bat_ind == batch_name[batches]))
      LM_Results <- lm(response~logrc+loglg+bgvector)
      Gamma[batches] <- LM_Results$coefficients[[1]]
      Alpha[batches] <- LM_Results$coefficients[[2]]
      Beta[batches] <- LM_Results$coefficients[[3]]
      Theta[batches] <- LM_Results$coefficients[[4]]
    }
    
    for(batches in 1:N_batch){
      Nu[batches] <- mean(bg_av[,batches])
      Sigma[batches] <- bgsq_av[batches] - Nu[batches]^2
    }
    Lambda <- sweep(sweep(Mu, 1, gene_len, FUN='*')%*%Group_Matrix, 2, tot_read, FUN = '*')
    for(batches in 1:N_batch){
      Pi[, bat_ind == batch_name[batches]] <- Gamma[batches] + sweep(sweep(Pi, 2, Alpha[batches]*log(tot_read[bat_ind==batch_name[batches]], FUN="+")), 1, Beta[batches]*log(gene_len), FUN="+")
    }
    cat("iteration", j)
    j <- j + 1
  }
  Z <- MH$Z/K
  Lambda_bg <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(batches in 1:N_batch){
    Lambda_bg[, bat_ind == batch_name[batches]] <- sweep(Lambda[, bat_ind == batch_name[batches]], 1, exp(bg_av[,batches]), FUN='*')
  }
  
  meanNu <- mean(Nu)
  bg_av = bg_av - meanNu
  Nu = Nu - meanNu
  Mu = Mu * exp(meanNu)
  for(batches in 1:N_batch){
    Gamma[batches] = Gamma[batches] + Theta[batches]*meanNu
  }
  
  # dropout imputation
  cat("start imputing ... ...")
  Y_impute <- matrix(0, ncol=Cell_N, nrow=Gene_N)
  for(cells in 1:Cell_N){
    for(genes in 1:Gene_N){
      Y_cg <- Y[genes, cells]
      lambda_bcg <- Lambda_bg[genes, cells]
      lambda_cg <- Lambda[genes, cells]
      lowerp <- ppois(q=Y_cg - 1, lambda=lambda_bcg)
      upperp <- ppois(q=Y_cg, lambda=lambda_bcg)
      lowerq <- qpois(p=lowerp, lambda=lambda_cg)
      upperq <- qpois(p=upperp, lambda=lambda_cg)
      sumweight <- upperp - lowerp
      tempsum <- 0
      for(imp in lowerq:upperq){
        low1p <- max(ppois(q=imp-1, lambda=lambda_cg), lowerp)
        up1p <- min(ppois(q=imp, lambda=lambda_cg), upperp)
        tempsum <- tempsum + imp*(up1p - low1p)
      }
      Y_impute[genes, cells] = tempsum/sumweight
    }
  }
  Y_impute <- Y_impute * Z + Lambda * (1 - Z)
  return(list(Mu=Mu, Alpha=Alpha, Beta=Beta, Gamma=Gamma, Theta=Theta, Nu=Nu, Sigma=Sigma, 
              Y_impute = Y_impute, bg=bg_av, Groups=group_name, Batches=batch_name))
}