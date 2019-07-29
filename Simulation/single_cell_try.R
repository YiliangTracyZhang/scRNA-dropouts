fitDABB <- function(Y, tot_read, bat_ind, bio_ind, gene_len, tech_para){
        G_num <- length(Y[ ,1])
        C_num <- length(Y[1, ])
        # former part seems to overlap with singlecell.R
        bat_num <- max(bat_ind)
        bio_num <- max(bio_ind)
        
        Y_binary <- ifelse(Y > 0, 1, 0)
        bat_mat <- c()
        for (i in 1:bat_num){
                bat_mat <- cbind(bat_mat, ifelse(bat_ind == i, 1, 0))
        }
        bat_cellnum <- colSums(bat_mat)
        bio_mat <- c()
        for (i in 1:bio_num){
                bio_mat <- cbind(bio_mat, ifelse(bio_ind == i, 1, 0))
        }
        
        RL_mat <- gene_len %*% t(tot_read) / 10^9
       # RL_mat <- gene_len %*% t(tot_read) 
        # adjustment of RL_mat/tot_read/gene_len may not be necessary in simulation data
        
        lnR_vec <- c()
        for (i in 1:C_num){
                 lnR_vec <- c(lnR_vec, rep(log(tot_read[i] / 10^6), G_num))
                #lnR_vec <- c(lnR_vec, rep(log(tot_read[i] ), G_num))
        }
         lnL_vec <- rep(log(gene_len / 10^3), C_num)
        # lnL_vec <- rep(log(gene_len ), C_num)
        
        Design_mat <- cbind(rep(1, G_num * C_num), lnR_vec, lnL_vec) # G*C,3
        Design_cov <- t(Design_mat) %*% Design_mat # 3*3
        
        iter_num <- tech_para$iternum
        max_error <- tech_para$error
        MH_num <- tech_para$mhnum
        sig_jump <- tech_para$jump
        dec_rate <- tech_para$jump_rate
        burn_num <- tech_para$burnin
        sample_num <- MH_num - burn_num
        
        iter <- 0
        loglike_old <- 1
        nu_Bat_old <- 0
        error <- 1
        ##sigma2_Bat_old <- 1
        
        #initial value
        nu_Bat <- rep(0, bat_num)
        sigma2_Bat <- rep(1, bat_num)
        coef <- c(0, 0.1, 0.1, 1)
        mu_mat <- mu_update(Y, bio_mat, RL_mat, matrix(1, G_num, C_num), rep(1, C_num))
        
        while (iter < iter_num & error > max_error){
                iter <- iter + 1
                #E-step
                if (iter > 1 & iter < 3){
                        sig_jump <- dec_rate * sig_jump
                }
                
                b_sample_mat <- c()
                loglike_mat <- c()
                expect_eta <- 0
                expect_eta_b <- 0
                p_weight <- 0
                Phi_mat <- 0
                ##loglike_all <- 0
                jump_mat <- c()
                
                b_old <- rnorm(C_num, as.vector(bat_mat %*% nu_Bat), 
                               as.vector(bat_mat %*% sigma2_Bat^0.5))
                lst_post_old <- b_post_loglike(Y, b_old, mu_mat, bio_mat, RL_mat, Y_binary,
                                               coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat)
                #loglike1 = loglike_cell + b_loglike_vec, loglike2 = loglike_cell, 
                #etaexpect = expect_eta, pmat = p_weight, Phimat = p_nodrop_mat,
                #loglikemat = loglike_mat
                log_post_old <- lst_post_old$loglike1
                
                log_post_bii <- lst_post_old$loglike2
                p_weightii <- lst_post_old$pmat
                expect_etaii <- lst_post_old$etaexpect
                Phi_matii <- lst_post_old$Phimat
                ##loglike_allii <- lst_post_old$loglikemat
                
                for (ii in 1:MH_num){
                        #Metropolis-Hasting
                        b_new <- rnorm(C_num, b_old, rep(sig_jump, C_num))
                        lst_post_new <- b_post_loglike(Y, b_new, mu_mat, bio_mat, RL_mat, Y_binary,
                                                       coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat)
                        log_post_new <- lst_post_new$loglike1
                        log_post_diff <- log_post_new - log_post_old
                        r <- ifelse(log_post_diff >= 0, 1, exp(log_post_diff))
                        unif <- runif(C_num)
                        b_old <- ifelse(r > unif, b_new, b_old)
                        log_post_old <- ifelse(r > unif, log_post_new, log_post_old)
                        jump_vec <- as.vector(ifelse(r > unif, 1, 0)) # what this means
                        jump_mat <- cbind(jump_mat, jump_vec)
                        Phi_matii <- t((1 - jump_vec) * t(Phi_matii) + jump_vec * t(lst_post_new$Phimat))
                        log_post_bii <- (1 - jump_vec) * log_post_bii + jump_vec * lst_post_new$loglike2
                        p_weightii <- t((1 - jump_vec) * t(p_weightii) + jump_vec * t(lst_post_new$pmat))
                        expect_etaii <- t((1 - jump_vec) * t(expect_etaii) +
                                                  jump_vec * t(lst_post_new$etaexpect))
                        ##loglike_allii <- t((1 - jump_vec) * t(loglike_allii) +
                        ##                   jump_vec * t(lst_post_new$loglikemat))
                        
                        #Monte Carlo
                        if (ii > burn_num){
                                # not mean()?
                                b_sample_mat <- cbind(b_sample_mat, b_old)
                                expect_eta <- expect_eta + expect_etaii
                                expect_eta_b <- expect_eta_b + t(t(expect_etaii) * as.vector(b_old))
                                p_weight <- p_weight + p_weightii
                                loglike_mat <- cbind(loglike_mat, log_post_bii)
                                Phi_mat <- Phi_mat + Phi_matii
                                ##loglike_all <- loglike_all + loglike_allii
                        }
                }
                expect_b <- rowMeans(b_sample_mat)
                expect_expb <- rowMeans(exp(b_sample_mat))
                expect_b2 <- rowMeans(b_sample_mat^2)
                expect_eta <- expect_eta / sample_num
                expect_eta_b <- expect_eta_b / sample_num
                p_weight <- p_weight / sample_num
                loglike_vec <- rowMeans(loglike_mat)
                loglike_new <- sum(loglike_vec)
                Phi_mat <- Phi_mat / sample_num
                ##loglike_all <- loglike_all / sample_num
                
                print('mean_jump_mat',mean(jump_mat))
                
                #M-step
                mu_mat <- mu_update(Y, bio_mat, RL_mat, p_weight, expect_expb)
                A1_12 <- t(Design_mat) %*% as.vector(rep(1, G_num) %*% t(expect_b))
                A1_22 <- G_num * sum(expect_b2)
                A1 <- cbind(rbind(Design_cov, t(A1_12)), rbind(A1_12, A1_22))
                A2 <- rbind(t(Design_mat) %*% as.matrix(as.vector(expect_eta)), 
                            sum(expect_eta_b)) 
                coef <- solve(A1) %*% A2
                expect_b_shift <- expect_b - mean(expect_b)
                
                nu_Bat <- (t(bat_mat) %*% expect_b_shift) / bat_cellnum
                sigma2_Bat <- (t(bat_mat) %*% ((expect_b - bat_mat %*% nu_Bat)^2)) / bat_cellnum
                
                print(c('iteration',iter,'new_log_like', loglike_new))
                
                error1 <- abs((loglike_new - loglike_old) / (loglike_old + 1))
                loglike_old <- loglike_new
                error2 <- 1e3 * mean((nu_Bat - nu_Bat_old)^2 / (nu_Bat_old + 1))
                nu_Bat_old <- nu_Bat
                error <- min(error1, error2)
                ##print(c(error1, error2))
        }
        ##loglike_gene <- rowSums(loglike_all)
        return(list(nu = nu_Bat, sigma2 = sigma2_Bat, mu = mu_mat, bsample = b_sample_mat, 
                    pweight = p_weight, coef = coef, eta = expect_eta, Phi = Phi_mat))
}

#update functions

mu_update <- function(Y, bio_mat, RL_mat, p_weight, exp_b){
        return((p_weight * Y) %*% bio_mat / (t(t(p_weight * RL_mat) * exp_b) %*% bio_mat))
}

b_post_loglike <- function(Y, b_vec, mu_mat, bio_mat, RL_mat, Y_binary,
                           coef, Design_mat, bat_mat, nu_Bat, sigma2_Bat){
        G_num <- length(Y_binary[,1])
        C_num <- length(Y_binary[1,])
        lambda_mat <- t(t(mu_mat %*% t(bio_mat) * RL_mat) * exp(b_vec))
        eta_mean_vec <- as.vector(cbind(Design_mat, b_vec) %*% coef)
        eta_mean <- matrix(eta_mean_vec, G_num, C_num)
        p_nodrop_mat <- matrix(pnorm(eta_mean_vec), G_num, C_num)
        zero_poisson <- p_nodrop_mat * (1 - Y_binary) * exp(-lambda_mat)
        zero_prob <- zero_poisson + (1 - p_nodrop_mat) * (1 - Y_binary)
        loglike_mat <- Y * log((1e-10) + lambda_mat) - lambda_mat * Y_binary + 
                log(p_nodrop_mat * Y_binary + zero_prob) 
        loglike_cell <- colSums(loglike_mat)
        b_loglike_vec <- - 1 / (2 * bat_mat %*% sigma2_Bat) * (b_vec - bat_mat %*% nu_Bat)^2
        p_weight <- Y_binary + zero_poisson / (zero_prob + Y_binary) 
        phi_mat <- dnorm(eta_mean) 
        expect_eta <- eta_mean + p_weight * phi_mat / p_nodrop_mat - 
                (1 - p_weight) * phi_mat / (1 - p_nodrop_mat)
        return(list(loglike1 = loglike_cell + b_loglike_vec, loglike2 = loglike_cell, 
                    etaexpect = expect_eta, pmat = p_weight, Phimat = p_nodrop_mat,
                    loglikemat = loglike_mat))
}

#####################################
##read data, we use alpha=0.1 here##
#####################################

raw.data <- read.table('/home/kl764/project/singlecell/simulation/alterA/0.1/read1.txt', sep = '', header = F)
data_all <- as.matrix(raw.data)
bio.group <- c(rep(1,2000))
batch.info <- c(rep(1,1000),rep(2,1000))
#total.count <- apply(raw.data[,1:2000],2,sum)
count.table <- read.table('/home/kl764/project/singlecell/simulation/alterA/0.1/R1.txt', sep = '', header = F)
total.count <- as.numeric(unlist(count.table))
gene_length <- as.numeric(as.vector(data_all[,2001]))
Y <- as.matrix(data_all[,c(1:2000)])
G_num <- length(Y[,1])
C_num <- length(Y[1,])
Y <- matrix(as.numeric(Y), G_num, C_num)

tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                  jump = 0.03, burnin = 50)
# source('/home/kl764/project/singlecell/git/scRNA-dropouts/fitDABB_function.R')

results <- fitDABB(Y, total.count, batch.info, bio.group, gene_length, tech_para)
write.table( results, paste0(('/home/kl764/project/singlecell/simulation/alterA/try-result.txt'), quote = F))
             