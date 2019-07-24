## Repeat 100 times for each setting
rm(list=ls())
set.seed(19720)

# the number of genes
p = 10000 
# the number of cells, each batch contains half of the cells
n = 2000 
#since genes we detect are fixed, so length of each gene is unchangable 
l = exp(rnorm(p))*100 

################################
## Scenario 1: One cell type ##
################################

# true expression level in each gene is fixed(under one cell type) 
mu = (rnorm(p, sd=0.5))^2 

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


## Altering Beta, when alpha=0.5, gamma=0, nu1=-1, nu2=1
A = Alpha[3]
C = Gamma[3]
Nu1 = nu1[3]
Nu2 = nu2[3]

for(j in 1:5){
        B = Beta[j]
        for(i in 1:100){
                # 1.modeling of non-dropout
                r = exp(rnorm(n)) # total read counts in each cell
                b1 = rnorm(1, mean = Nu1, sd = Sigma) # batch effect in batch 1
                b2 = rnorm(1, mean = Nu2, sd = Sigma) # batch effect in batch 2
                b = c(rep(b1, n/2), rep(b2, n/2))
                exp.b = c(rep(exp(b1), n/2), rep(exp(b2), n/2)) # exponential batch effect in each cell
                #  Lambda=(mu*l*b) %*% t(l) 
                Lambda = (mu*l) %*% t(r*exp.b) # non-dropout read count
                
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
      #  write.table(Y, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterB/B', B, '/read', j, '.txt'), quote = F, col.names = F, row.names = F)
        #file_name <- paste0('~/Desktop/simulation-data',B,'txt')
        # write.table(Y[[B]],file_name)
        ## Repeat 100 times for each setting
        write.table(Y, paste0('/home/kl764/project/singlecell/simulation/alterB/', B, '/read', i,'.txt'), 
                    quote = F, col.names = F, row.names = F)
        }
}