## Repeat 100 times for each setting
rm(list=ls())

##################################
## Scenario 2: Three cell types ##
##################################

# the number of genes
p = 10000 
# the number of cells, each batch contains half of the cells
n = 2000 
# parameters
#Alpha = c(0,0.2,0.5,0.7,1)
Alpha = c(-0.1,0,0.1,0.3,0.5)
#Beta = c(0,0.2,0.5,0.7,1)
Beta = c(-0.1,0,0.1,0.3,0.5)
Gamma = c(-2,-1,0,1,2)
nu1 = c(-2,-1.5,-1,-0.5,0)
nu2 = c(2,1.5,1,0.5,0)
Sigma = 1
Theta = 1

## Altering Nu1, when alpha=0.5, beta=0.5, gamma=0
A = Alpha[3]
B = Beta[3]
C = Gamma[3]

for(j in 1:5){
        Nu1 = nu1[j]  
        Nu2 = nu2[j]
        for(i in 1:100){
                # 1.modeling of non-dropout
                mu1 = exp(rnorm(p))*5 # true expression level in each gene(cell type 1)
                mu2 = exp(rnorm(p))*5 # true expression level in each gene(cell type 2)
                mu3 = exp(rnorm(p))*5 # true expression level in each gene(cell type 3)
                r = exp(rnorm(n))*5 # total read counts in each cell
                l = exp(rnorm(p))*5 # length of each gene
                b1 = rnorm(1, mean = Nu1, sd = Sigma) # batch effect in batch 1
                b2 = rnorm(1, mean = Nu2, sd = Sigma) # batch effect in batch 2
                b = c(rep(exp(b1), n/2), rep(exp(b2), n/2)) # batch effect in each cell
                #  Lambda=(mu*l*b) %*% t(l) 
                #  Lambda = (mu*l) %*% t(r*b) # non-dropout read count
                lam1 = (mu1*l) %*% t((r*b)[1:250]) 
                lam2 = (mu2*l) %*% t((r*b)[251:600]) 
                lam3 = (mu3*l) %*% t((r*b)[601:1000]) 
                lam4 = (mu1*l) %*% t((r*b)[1001:1300]) 
                lam5 = (mu2*l) %*% t((r*b)[1301:1650]) 
                lam6 = (mu3*l) %*% t((r*b)[1651:2000]) 
                Lambda <- cbind(lam1,lam2,lam3,lam4,lam5,lam6)
                # 2.modeling of dropout effect
                # Epsilon <- matrix(rnorm(n*p), n, p)
                # Z=Epsilon<=(C+(A*log(r)+b) %*% t(rep(1,p))+B*rep(1,n) %*% t(l))
                Epsilon <- matrix(rnorm(p*n), p, n)
                Z = Epsilon <= (
                        C*rep(1,p)%*%t(rep(1,n))
                        + rep(1,p) %*% t(A*log(r)+b)
                        + log(l) %*% t(B*rep(1,n))
                ) 
                # 3.expression of read count
                #Y=Z*matrix(rpois(n*p, as.vector(Lambda)), n, p)
                Y = Z*matrix(rpois(p*n, as.vector(Lambda)), p, n)
        }
        
        #  write.table(Y, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterNu/Nu', Nu1, '/read', j, '.txt'), 
        #               quote = F, col.names = F, row.names = F)
        write.table(Y, paste0('/home/kl764/project/singlecell/simulation/alterNu.ct/', j, '.txt'), 
                    quote = F, col.names = F, row.names = F)
}
