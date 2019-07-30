## Repeat 100 times for each setting
rm(list=ls())
set.seed(19720)

# the number of genes
p = 10000
# the number of cells, each batch contains half of the cells
n = 1000
#since genes we detect are fixed, so length of each gene is unchangable 
l = exp(rnorm(p))*100 
write.table(l, paste0('/home/kl764/project/singlecell/simulation/alterB/l.txt'), 
            quote = F, col.names = F, row.names = F)

################################
## Scenario 1: One cell type ##
################################

# true expression level in each gene is fixed(under one cell type) 
mu = (rnorm(p, sd=0.5))^2 
write.table(mu, paste0('/home/kl764/project/singlecell/simulation/alterB/mu.txt'), 
            quote = F, col.names = F, row.names = F)

# parameters
Alpha = c(0,0.05,0.10,0.15,0.20)
Beta = c(0,0.05,0.10,0.15,0.20)
Gamma = c(-1.5,-1,-0.5,0,0.5)
nu1 = c(-0.2,-0.15,-0.10,-0.05,0)
nu2 = c(0.2,0.15,0.10,0.05,0)
Sigma = 1
Theta = 1


## Altering Beta, when alpha=0.1, gamma=-0.5, nu1=-0.1, nu2=0.1
A = Alpha[3]
C = Gamma[3]
Nu1 = nu1[3]
Nu2 = nu2[3]

for(j in 1:5){
        B = Beta[j]
        for(i in 1:100){
                # 1.modeling of non-dropout
                r = exp(rnorm(n)) # total read counts in each cell
                b1 = rnorm(n/2, mean = Nu1, sd = Sigma) # batch effect in batch 1
                b2 = rnorm(n/2, mean = Nu2, sd = Sigma) # batch effect in batch 2
                b = c(b1, b2)
                exp.b = c(exp(b1), exp(b2))  # exponential batch effect in each cell
                Lambda = (mu*l) %*% t(r*exp.b) # non-dropout read count
                
                # 2.modeling of dropout effect
                Epsilon <- matrix(rnorm(p*n), p, n)
                Z = Epsilon <= (
                        C*rep(1,p)%*%t(rep(1,n))
                        + rep(1,p) %*% t(A*log(r)+Theta*b)
                        + log(l) %*% t(B*rep(1,n))
                ) 
                
                # 3.expression of read count
                Y = cbind(Z*matrix(rpois(p*n, as.vector(Lambda)), p, n),l)
                
        write.table(Y, paste0('/home/kl764/project/singlecell/simulation/alterB/', B, '/read', i,'.txt'), 
                   quote = F, col.names = F, row.names = F)
        write.table(r, paste0('/home/kl764/project/singlecell/simulation/alterB/', B, '/R', i,'.txt'), 
                    quote = F, col.names = F, row.names = F)
        }
}