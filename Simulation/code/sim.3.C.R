## Repeat 100 times for each setting
rm(list=ls())
set.seed(19720)

# the number of genes
p = 10000 
# the number of cells, each batch contains half of the cells
n = 2000 
#since genes we detect are fixed, so length of each gene is unchangable 
l = exp(rnorm(p))*100 

# parameters
#Alpha = c(0,0.2,0.5,0.7,1)
Alpha = c(0,0.05,0.10,0.15,0.20)
#Beta = c(0,0.2,0.5,0.7,1)
Beta = c(0,0.05,0.10,0.15,0.20)
# Gamma = c(-2,-1,0,1,2)
Gamma = c(-1.5,-1,-0.5,0,0.5)
#nu1 = c(-2,-1.5,-1.0,-0.5,0)
nu1 = c(-0.2,-0.15,-0.10,-0.05,0)
#nu2 = c(2,1.5,1.0,0.5,0)
nu2 = c(0.2,0.15,0.10,0.05,0)
Sigma = 1
Theta = 1

##################################
## Scenario 2: Three cell types ##
##################################

# note that each cell type has its own expression pattern
# assuming in batch 1: type 1 has 250 cells, type 2 has 350 cells, type 3 has 400 cells
# in batch 2:  type 1 has 300 cells, type 2 has 350 cells, type 3 has 350 cells
mu1 = (rnorm(p))^2 # true expression level in each gene(cell type 1)
mu2 = (rnorm(p))^2 # true expression level in each gene(cell type 2)
mu3 = (rnorm(p))^2 # true expression level in each gene(cell type 3)

## Altering Gamma, when alpha=0.5, beta=0.5, nu1=-1, nu2=-1
A = Alpha[3]
B = Beta[3]
Nu1 = nu1[3]
Nu2 = nu2[3]

for(j in 1:5){
        C = Gamma[j]
        for(i in 1:100){
                # 1.modeling of non-dropout
                r = exp(rnorm(n)) # total read counts in each cell
                b1 = rnorm(1, mean = Nu1, sd = Sigma) # batch effect in batch 1
                b2 = rnorm(1, mean = Nu2, sd = Sigma) # batch effect in batch 2
                b = c(rep(b1, n/2), rep(b2, n/2))
                exp.b = c(rep(exp(b1), n/2), rep(exp(b2), n/2)) # exponential batch effect in each cell
                #  Lambda=(mu*l*b) %*% t(l) 
                # Lambda = (mu*l) %*% t(r*b) # non-dropout read count
                lam1 = (mu1*l) %*% t((r*exp.b)[1:250]) 
                lam2 = (mu2*l) %*% t((r*exp.b)[251:600]) 
                lam3 = (mu3*l) %*% t((r*exp.b)[601:1000]) 
                lam4 = (mu1*l) %*% t((r*exp.b)[1001:1300]) 
                lam5 = (mu2*l) %*% t((r*exp.b)[1301:1650]) 
                lam6 = (mu3*l) %*% t((r*exp.b)[1651:2000]) 
                Lambda <- cbind(lam1,lam2,lam3,lam4,lam5,lam6)
                
                # 2.modeling of dropout effect
                # Epsilon <- matrix(rnorm(n*p), n, p)
                # Z=Epsilon<=(C+(A*log(r)+b) %*% t(rep(1,p))+B*rep(1,n) %*% t(l))
                Epsilon <- matrix(rnorm(p*n), p, n)
                Z = Epsilon <= (
                        C*rep(1,p)%*%t(rep(1,n))
                        + rep(1,p) %*% t(A*log(r)+b*Theta)
                        + log(l) %*% t(B*rep(1,n))
                ) 
                
                # 3.expression of read count
                #Y=Z*matrix(rpois(n*p, as.vector(Lambda)), n, p)
                Y = Z*matrix(rpois(p*n, as.vector(Lambda)), p, n)
                Y = cbind(Y,l)
                
                # write.table(Y, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterC/C', C, '/read', j, '.txt'), 
                #             quote = F, col.names = F, row.names = F)
                write.table(Y, paste0('/home/kl764/project/singlecell/simulation/alterC.ct/',  C, '/read', i, '.txt'), 
                            quote = F, col.names = F, row.names = F)
        }

}
