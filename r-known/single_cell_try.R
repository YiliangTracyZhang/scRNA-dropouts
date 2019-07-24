#Read the data, we use deng's data as a test of our method.
#data_all <- as.matrix(new.processed.data)
#Scenario = args[1]
#raw.data <- read.table('/home/kl764/project/singlecell/', sep = '', header = F)
#raw.data <- read.table('~/Documents/singlecell/simulation/extreme.txt', sep = '', header = F)
#setwd('~/Documents/singlecell/simulation/alterA-r-known/0.2/')
setwd('/home/kl764/project/singlecell/git/scRNA-dropouts/r-known/alterA-r-known/0.2/')
#data <- dir(path = '~/Documents/singlecell/simulation/alterA-r-known/0.2/', pattern = ".txt")
data <- dir(path = '/home/kl764/project/singlecell/git/scRNA-dropouts/r-known/alterA-r-known/0.2/', pattern = ".txt")
parameters<-c()
for( i in 1:length(data)){
raw.data<- read.table(data[i], sep = '', header = F)
data_all <- as.matrix(raw.data)
#head(data_all)
#bio.group <- c(rep(1,35), rep(2,19))
#bio.group <- c(rep(1,250), rep(2,350),rep(3,400),rep(1,300),rep(2,350),rep(3,350))
bio.group <- c(rep(1,100))
#bio.group <- c(rep(1,1000))
#batch.info <- c(rep(1,27),rep(2,8),rep(1,10),rep(2,9))
batch.info <- c(rep(1,50),rep(2,50))
#batch.info <- c(rep(1,500),rep(2,500))
#total.count <- apply(new.processed.data[,3:56],2,sum)
total.count <- apply(raw.data[,1:100],2,sum)
#total.count <- apply(raw.data[,1:1000],2,sum)
#Y <- as.matrix(data_all[,c(3:56)])
Y <- as.matrix(data_all[,c(1:100)])
#Y <- as.matrix(data_all[,c(1:1000)])
G_num <- length(Y[,1])
C_num <- length(Y[1,])
Y <- matrix(as.numeric(Y), G_num, C_num)
gene_length <- as.numeric(as.vector(data_all[,101]))
#gene_length <- as.numeric(as.vector(data_all[,1001]))
tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                  jump = 0.03, burnin = 50)
source('/home/kl764/project/singlecell/git/scRNA-dropouts/fitDABB_function.R')
#source('/Users/kexuanliang/documents/singlecell/git/scRNA-dropouts/fitDABB_function.R')

#fit the model:
#Y is the G*C read count matrix, total.count is a vector of each cell's total read counts,
#batch.info is a vector, the index of the batch each cell belongs to, bio.group is a 
#vector, the index of the biological group each cell belongs to, gene_length is a vector,
#length of each RNA, and tech_para is a group of parameters that are used to fit the model.
#Among these parameters, jump_rate and jump are used to control the acceptance rate of 
#Metropolis-Hasting alogrithm. A recommanded rate ranges from 0.2 to 0.45.
results <- fitDABB(Y, total.count, batch.info, bio.group, gene_length, tech_para)
#write.table(results, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterA-r-known.result/0.2/',i,'.txt'), 
#            quote = F, col.names = F, row.names = F)
write.table(results, paste0('/home/kl764/project/singlecell/git/scRNA-dropouts/r-known/alterA-r-known/0.2-result/',i,'.txt'), 
            quote = F, col.names = F, row.names = F)
parameters <- cbind(parameters,rbind(results$nu,results$sigma2,results$coef))
}

para<-data.frame(t(as.matrix(parameters)))
names(para)<-c('nu','sigma2','gamma','alpha','beta','thita')
#write.table(para, paste0('/Users/kexuanliang/documents/singlecell/simulation/alterA-r-known.result/0.2.txt'))
write.table(para, paste0('/home/kl764/project/singlecell/git/scRNA-dropouts/r-known/alterA-r-known/0.2-result/para.txt'))

p<-boxplot(para)
pdf('/home/kl764/project/singlecell/git/scRNA-dropouts/r-known/alterA-r-known/0.2-result/boxplot.pdf')
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