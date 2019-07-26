## Description: this R script is written for testing simulation result
## we compare 2 scenarios here, each includes 100 iteration data, detailed info as described below
## simulation parameter setting:
### Alpha = c(0,0.05,0.10,0.15,0.20)
### Beta = c(0,0.05,0.10,0.15,0.20)
### Gamma = c(-1.5,-1,-0.5,0,0.5)
### nu1 = c(-0.2,-0.15,-0.10,-0.05,0)
### nu2 = c(0.2,0.15,0.10,0.05,0)
### theta = 1, sigma = 1 

## NOTE 1: 
## without the [simulated] information of total count, we checked the estimated parameters under different datasets,
## results were much poorer (not shown in this script)
## NOTE 2: 
## magnifying the parameters when simulating data DOES NOT improve the results ( with/ without R info), e.g. setting A = 1, B = 1, C = -2, Nu1 = -1,
## so larger parameters may not improve the properties of our model
## NOTE 3:
## input n=100, p=100, the time for running the model should be around 20 - 30 min under 1 cell type. Longer time taken (30-45min) under 3 cell types.
## but input n=2000, p=10000, the time become much longer even running on cluster


###############
##Scenario 1 ##
###############

## simulated dataset description:
## n = 100, p = 100
## cell type: 1
## batch: 2, each contains 50 cells
## Altering Alpha, when Beta=0.5, Gamma=0, Nu1=-1, Nu2=1
## inputs are two extreme cases of A, i.e. A = 0 and A = 0.2

# input simulation data of [total count of cells] 
setwd('~/Documents/singlecell/simulation/alterA/Rinfo-0.2/')
count <- dir(path = '~/Documents/singlecell/simulation/alterA/Rinfo-0.2/', pattern = ".txt")
for( i in 1:length(count)){
        count.table[i] <- read.table(count[i], sep ='', header = F, as.is = T)
}

# input simulation data of [read count matrix] 
setwd('~/Documents/singlecell/simulation/alterA/0.2/')
data <- dir(path = '~/Documents/singlecell/simulation/alterA/0.2/', pattern = ".txt")

parameters1 <-c()

# fitting data into our model 
for( i in 1:length(data)){
        
        total.count <- as.numeric(unlist(count.table[i]))
        #total.count <- apply(raw.data[,1:100],2,sum)
        
        raw.data <- read.table(data[i], sep = '', header = F)
        
        gene_length <- as.numeric(as.vector(data_all[,101]))   
        
        data_all <- as.matrix(raw.data)

        bio.group <- c(rep(1,100))

        batch.info <- c(rep(1,50),rep(2,50))
        
        Y <- as.matrix(data_all[,c(1:200)])
        G_num <- length(Y[,1])
        C_num <- length(Y[1,])
        Y <- matrix(as.numeric(Y), G_num, C_num)
       
        
        tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                          jump = 0.03, burnin = 50)
        
        #source('/home/kl764/project/singlecell/git/scRNA-dropouts/fitDABB_function.R')
        source('/Users/kexuanliang/documents/singlecell/git/scRNA-dropouts/fitDABB_function.R')
        
        #fit the model:
        #Y is the G*C read count matrix, total.count is a vector of each cell's total read counts,
        #batch.info is a vector, the index of the batch each cell belongs to, bio.group is a 
        #vector, the index of the biological group each cell belongs to, gene_length is a vector,
        #length of each RNA, and tech_para is a group of parameters that are used to fit the model.
        #Among these parameters, jump_rate and jump are used to control the acceptance rate of 
        #Metropolis-Hasting alogrithm. A recommanded rate ranges from 0.2 to 0.45.
        
        results <- fitDABB(Y, total.count, batch.info, bio.group, gene_length, tech_para)
        
        #write.table(results, paste0('/Users/kexuanliang/documents/singlecell/simulation/para-result/',i,'.txt'), 
        #   quote = F, col.names = F, row.names = F)
        
        parameters1 <- cbind(parameters,rbind(results$nu,results$sigma2,results$coef))
        
}

# store the parameter information into a data.frame
para1<-data.frame(t(as.matrix(parameters1)))
names(para1)<-c('nu1','nu2','sigsq1','sigsq2','gamma','alpha','beta','theta')
write.table(para1, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterA/0.2/para-results.txt'), 
   quote = F, row.names = F)

# visualize our result
p<-boxplot(para1)
png(file="boxplot")
dev.off()

###############
##Scenario 2 ##
###############
## simulated dataset description:
## n = 200, p = 100
## batch: 2, each contains 100 cells
## cell type: 3, batch 1 (25,35,40), batch 2 (30,35,35)
## Altering Alpha, when Beta=0.5, Gamma=0, Nu1=-1, Nu2=1
## input is an extreme case of A, A = 0.2


setwd('~/Documents/singlecell/simulation/alterA.ct/Rinfo-0.2/')
count <- dir(path = '~/Documents/singlecell/simulation/alterA.ct/Rinfo-0.2/', pattern = ".txt")
for( i in 1:length(count)){
        count.table[i] <- read.table(count[i], sep ='', header = F, as.is = T)
}

setwd('~/Documents/singlecell/simulation/alterA.ct/0.2/')
data <- dir(path = '~/Documents/singlecell/simulation/alterA.ct/0.2/', pattern = ".txt")

parameters2<-c()

for( i in 1:length(data)){
        
        raw.data <- read.table(data[i], sep = '', header = F)
        
        total.count <- as.numeric(unlist(count.table[i]))
        #total.count <- apply(raw.data[,1:200],2,sum)
        
        data_all <- as.matrix(raw.data)
        
        gene_length <- as.numeric(as.vector(data_all[,201]))
        
        bio.group <- c(rep(1,25), rep(2,35),rep(3,40),rep(1,30),rep(2,35),rep(3,35))
        
        batch.info <- c(rep(1,100),rep(2,100))
        
        Y <- as.matrix(data_all[,c(1:200)])
        G_num <- length(Y[,1])
        C_num <- length(Y[1,])
        Y <- matrix(as.numeric(Y), G_num, C_num)
        
        
        tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                          jump = 0.03, burnin = 50)
        
        #source('/home/kl764/project/singlecell/git/scRNA-dropouts/fitDABB_function.R')
        source('/Users/kexuanliang/documents/singlecell/git/scRNA-dropouts/fitDABB_function.R')
     
        results <- fitDABB(Y, total.count, batch.info, bio.group, gene_length, tech_para)
        
        #write.table(results, paste0('/Users/kexuanliang/documents/singlecell/simulation/para-result/',i,'.txt'), 
        #   quote = F, col.names = F, row.names = F)
        
        parameters2 <- cbind(parameters,rbind(results$nu,results$sigma2,results$coef))
}

para2<-data.frame(t(as.matrix(parameters2)))
names(para2)<-c('nu1','nu2','sigsq1','sigsq2','gamma','alpha','beta','theta')
write.table(para2, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterA.ct/0.2/para-results.txt'), 
            quote = F, row.names = F)

p<-boxplot(para2)
png(file="boxplot")
dev.off()


#library(ggplot2)
#a1.coef <- data.frame(1, 1:100, results$coef)
#Qualtiy Control, alternative hypothesis can be chosen as 'left', 'right' and 'two side'.
#DABB_QC(results, alternative = 'right')
#refit after deleting outliers
#results.qc <- fitDABB(Y[,-c(22, 33, 35)], total.count[-c(22, 33, 35)], batch.info[-c(22, 33, 35)], bio.group[-c(22, 33, 35)], 
 #                 gene_length, tech_para)
#Differential Expression
#pvl <- DABB_DE(Y[,-c(22, 33, 35)], total.count[-c(22, 33, 35)], bio.group[-c(22, 33, 35)], gene_len = gene_length, results.qc, 
 #                sample_num = 200)

#sort(pvl$p.value, index.return = T)
#sort(p.adjust(pvl$p.value))

#Visualization, method can be chosen as 'PCA' and 'ISOmap'
#library(ggplot2)
#library(vegan)
#visual <- DABB_visualize(Y[,-c(22,33,35)], total.count[-c(22,33,35)], gene_length, results$pweight,
#                        results$bsample, method = 'PCA', k = 10)
#visual <- DABB_visualize(Y[,-c(22,33,35)], total.count[-c(22,33,35)], gene_length, results.qc$pweight,
                     #   results.qc$bsample, method = 'PCA', k = 10)
##visual <- visual$points
#pc1 <- visual[,1]
#pc2 <- visual[,2]
#pc1.16.88 <- pc1[1:26]# cell 22 is deleted after quality control

#pc1.16.193 <- pc1[27:32]# cell 33 and 35 are deleted after quality control

#pc1.8.88 <- pc1[33:42]

#pc1.8.193 <- pc1[43:51]

#pc2.16.88 <- pc2[1:26]# cell 22 is deleted after quality control

#pc2.16.193 <- pc2[27:32]# cell 33 and 35 are deleted after quality control

#pc2.8.88 <- pc2[33:42]

#pc2.8.193 <- pc2[43:51]
        
#p <- ggplot()
#p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
#p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')
#p1 <- ggplot()
#p1 <- p1 + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
#p1 <- p1 + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
#p1 <- p1 + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
#p1 <- p1 + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
#p1 <- p1 + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')
####################################################################
#Visualization of the unfitted data.
#G_num <- length(Y[,1])
#C_num <- length(Y[1,])
#sample_num <- length(results$bsample[1,])
#visual <- DABB_visualize(Y, total.count, gene_length, matrix(1, G_num, C_num),
                       #  b_sample_mat = matrix(0, C_num, sample_num), method = 'PCA', k = 10)
##visual <- visual$points
#pc1 <- visual[,1]
#pc2 <- visual[,2]
#pc1.16.88 <- pc1[1:27]
#pc1.16.193 <- pc1[28:35]
#pc1.8.88 <- pc1[36:45]
#pc1.8.193 <- pc1[46:54]
#pc2.16.88 <- pc2[1:27]
#pc2.16.193 <- pc2[28:35]
#pc2.8.88 <- pc2[36:45]
#pc2.8.193 <- pc2[46:54]

#p <- ggplot()
#p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
#p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')


#Visualization of the unfitted data.
#G_num <- length(Y[,1])
#C_num <- length(Y[1,])
#sample_num <- length(results$bsample[1,])
#visual <- DABB_visualize(Y, total.count, gene_length, matrix(1, G_num, C_num),
                        # b_sample_mat = matrix(0, C_num, sample_num), method = 'PCA', k = 10)
##visual <- visual$points
#pc1 <- visual[,1]
#pc2 <- visual[,2]
#pc1.16.88 <- pc1[1:27]
#pc1.16.193 <- pc1[28:35]
#pc1.8.88 <- pc1[36:45]
#pc1.8.193 <- pc1[46:54]
#pc2.16.88 <- pc2[1:27]
#pc2.16.193 <- pc2[28:35]
#pc2.8.88 <- pc2[36:45]
#pc2.8.193 <- pc2[46:54]

#p <- ggplot()
#p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
#p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
#p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')