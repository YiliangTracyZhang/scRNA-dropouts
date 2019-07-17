## Repeat 100 times for each setting
## Scenario 1: One cell type
p=10000
n=2000

Alpha=c(0,0.2,0.5,0.7,1)
Beta=c(0,0.2,0.5,0.7,1)
Gamma=c(-2,-1,0,1,2)
nu1=c(-2,-1.5,-1,-0.5,0)
nu2=c(2,1.5,1,0.5,0)
Sigma=1
Theta=1
## alter alpha, beta=0.5, gamma=0, nu1=-1
B=Beta[3]
C=Gamma[3]
Nu1=nu1[3]
Nu2=nu2[3]
for(j in 1:5){
  A=Alpha[j]
  for(i in 1:100){
  # non-dropout
    mu=exp(rnorm(p))
    r=exp(rnorm(n))
    l=exp(rnorm(p))
    b1=rnorm(1, mean=Nu1, sd=Sigma)
    b2=rnorm(1, mean=Nu2, sd=Sigma)
    b=c(rep(exp(b1), n/2), rep(exp(b2), n/2))
    Lambda=(mu*l*b) %*% t(l) 

    # dropout
    Epsilon <- matrix(rnorm(n*p), n, p)
    Z=Epsilon<=(C+(A*log(r)+b) %*% t(rep(1,p))+B*rep(1,n) %*% t(l))
      
    # Expression
    Y=Z*matrix(rpois(n*p, as.vector(Lambda)), n, p)
  }
  write.table(Y, paste0('somewhere/alterA/A', A, '/read', j, '.txt'), quote = F, col.names = F, row.names = F)
}
## alter beta, alpha=0.5, gamma=0, nu1=-1
# non-dropout



# dropout

## alter gamma, alpha=0.5, beta=0.5, nu1=-1
# non-dropout



# dropout

## alter nu1, alpha=0.5, beta=0.5, gamma=0
# non-dropout



# dropout

## Scenario 2: Three cell types
# non-dropout



# dropout

