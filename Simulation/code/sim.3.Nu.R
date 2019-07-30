## Repeat 100 times for each setting
rm(list=ls())
set.seed(19720)

# the number of genes
p = 10000
# the number of cells, each batch contains half of the cells
n = 1000
#since genes we detect are fixed, so length of each gene is unchangable 
l = exp(rnorm(p))*100 
write.table(l, paste0('/home/kl764/project/singlecell/simulation/alterNu.ct/l.txt'), 
            quote = F, col.names = F, row.names = F)

##################################
## Scenario 2: Three cell types ##
##################################

# note that each cell type has its own expression pattern
# assuming in batch 1: type 1 has 250 cells, type 2 has 350 cells, type 3 has 400 cells
# in batch 2:  type 1 has 300 cells, type 2 has 350 cells, type 3 has 350 cells
mu1 = (rnorm(p))^2 # true expression level in each gene(cell type 1)
mu2 = (rnorm(p))^2 # true expression level in each gene(cell type 2)
mu3 = (rnorm(p))^2 # true expression level in each gene(cell type 3)
mu = rbind(cbind('type-1','type-2','type-3') ,cbind(mu1,mu2,mu3))
write.table(mu, paste0('/home/kl764/project/singlecell/simulation/alterNu.ct/mu.txt'), 
            quote = F, col.names = F, row.names = F)

## Altering Nu1, when alpha=0.1, beta=0.1, gamma=-0.5
A = Alpha[3]
B = Beta[3]
C = Gamma[3]

for(j in 1:5){
        Nu1 = nu1[j]  
        Nu2 = nu2[j]
        for(i in 1:100){
                
                # 1.modeling of non-dropout
                
                r = exp(rnorm(n))# total read counts in each cell
                b1 = rnorm(n/2, mean = Nu1, sd = Sigma) # batch effect in batch 1
                b2 = rnorm(n/2, mean = Nu2, sd = Sigma) # batch effect in batch 2
                b = c(b1, b2)
                exp.b = c(exp(b1), exp(b2))  # exponential batch effect in each cell
                # non-dropout read count
                lam1 = (mu1*l) %*% t((r*exp.b)[1:150]) 
                lam2 = (mu2*l) %*% t((r*exp.b)[151:300]) 
                lam3 = (mu3*l) %*% t((r*exp.b)[301:500]) 
                lam4 = (mu1*l) %*% t((r*exp.b)[501:700]) 
                lam5 = (mu2*l) %*% t((r*exp.b)[701:850]) 
                lam6 = (mu3*l) %*% t((r*exp.b)[851:1000]) 
                Lambda <- cbind(lam1,lam2,lam3,lam4,lam5,lam6)
                
                # 2.modeling of dropout effect
                
                Epsilon <- matrix(rnorm(p*n), p, n)
                Z = Epsilon <= (
                        C*rep(1,p)%*%t(rep(1,n))
                        + rep(1,p) %*% t(A*log(r)+Theta*b)
                        + log(l) %*% t(B*rep(1,n))
                ) 
                
                # 3.expression of read count
                
                Y = cbind(Z*matrix(rpois(p*n, as.vector(Lambda)), p, n),l)
                
                write.table(Y, paste0('/home/kl764/project/singlecell/simulation/alterNu.ct/',  Nu1, '/read', i, '.txt'), 
                           quote = F, col.names = F, row.names = F)
                write.table(r, paste0('/home/kl764/project/singlecell/simulation/alterNu.ct/',  Nu1, '/R', i, '.txt'), 
                            quote = F, col.names = F, row.names = F)
        }
}
