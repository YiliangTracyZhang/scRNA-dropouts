library(Biobase)
#library(biomaRt)
# mart = useMart('ensembl')
# listDatasets(mart)
# mart.gene = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#getENSG <- function(gene, mart = mart.gene) {
#  results <- getBM(attributes = c("ensembl_gene_id","start_position","end_position","transcript_length","hgnc_symbol"),
 #                  filters    = "hgnc_symbol", values = gene, mart = mart)
#  return(results$transcript_length)
#}


# read the expression set of usoskin data(filtered, we use zinger processed data here)
load("~/Documents/singlecell/git/scRNA-dropouts/real-dataset/esetUsoskin.RData")
usoskin.count <- eset
# gene_len <- getENSG(feature, mart.gene)
genes_len=read.table("/Users/kexuanliang/Documents/singlecell/git/scRNA-dropouts/real-dataset/usoskin/mm10_ccds_length.txt",stringsAsFactors=F)
colnames(genes_len) <- c("gene_name","gene_len")

pheno <- pData(usoskin.count)
bat_ind <- as.numeric(as.factor(pheno$`Picking sessions`)) # batches are indexed with 1(cold), 2(RT-1), 3(RT-2)
bio_ind <- pheno$`Level 1` # cell types are indexed with 1(NF), 2(NP), 3(PHP), 4(TH)
bio_ind <- gsub('NF', '1', bio_ind)
bio_ind <- gsub('NP', '2', bio_ind)
bio_ind <- gsub('PEP', '3', bio_ind)
bio_ind <- gsub('TH', '4', bio_ind)
bio_ind <- as.numeric(as.factor(bio_ind))
feature <- unlist(fData(usoskin.count))

expr.mat <- exprs(usoskin.count)
rownames(expr.mat) <- feature
expr.mat <- expr.mat[rownames(expr.mat) %in% genes_len$gene_name, ]
gene_len <- genes_len[ match(rownames(expr.mat), genes_len$gene_name), 2]

# tot_read <- unlist(pheno$Reads)
tot_read <- colSums(expr.mat)

C_num <- length(expr.mat[,1])
G_num <- length(expr.mat[1,])
Y <- expr.mat

tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                  jump = 0.03, burnin = 50)
source('/Users/kexuanliang/Documents/singlecell/git/scRNA-dropouts/fitDABB_function.R')

results <- fitDABB(Y, tot_read , bat_ind, bio_ind, gene_len, tech_para)
results$coef
results.para <- data.frame(c(as.numeric(results$nu), as.numeric(results$sigma2), as.numeric(results$coef)))
row.names(results.para) <- c('nu1','nu2', 'nu3', 'sig2-1', 'sig2-2', 'sig2-3', 'coef-1', 'coef-2', 'coef-3', 'coef-4' )
# df <- data.frame(unlist(results.list))
write.table(results.para, paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/results-2.txt'), quote = F, col.names = F)
write.table(data.frame(results$mu), paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/true.expression-2.txt'), quote = F)
write.table(data.frame(results$bsample), paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/bsample-2.txt'), quote = F, col.names = F)
write.table(data.frame(results$pweight), paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/pweight-2.txt'), quote = F, col.names = F)
write.table(data.frame(results$eta), paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/eta-2.txt'), quote = F, col.names = F)
write.table(data.frame(results$Phi), paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/phi-2.txt'), quote = F, col.names = F)


quality.control <- DABB_QC(results, alternative = 'right')
pvalue_out <- data.frame(cbind( 1:length(quality.control), quality.control))
names(pvalue_out) <- c('gene.num', 'p_value')
outliner <- data.frame(pvalue_out[pvalue_out$p_value <= 0.05,])
write.table(outliner,paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/qc.txt'), sep = '')

#Differential Expression without qc
pvl <- DABB_DE(Y, tot_read, bio_ind, gene_len, results, sample_num = 200) #not use
# refit after deleting outliers
#results <- fitDABB(Y[,-c(22, 33, 35)], total.count[-c(22, 33, 35)], batch.info[-c(22, 33, 35)], bio.group[-c(22, 33, 35)],
#                   gene_length, tech_para)

results.qc <- fitDABB(Y[,-outliner$gene.num], tot_read[-outliner$gene.num], bat_ind[-outliner$gene.num], bio_ind[-outliner$gene.num],
                      gene_len, tech_para)

#Differential Expression
pvl.qc <- DABB_DE(Y[-outliner$gene.num,], tot_read[-outliner$gene.num], bio_ind[-outliner$gene.num], gene_len, results.qc, sample_num = 200)#not use

#sort(pvl$p.value, index.return = T)
#sort(p.adjust(pvl$p.value))

#Visualization, method can be chosen as 'PCA' and 'ISOmap'
library(ggplot2)
library(vegan)

## visualize our result before quality control
visual <- DABB_visualize(Y, tot_read, gene_len, results$pweight,
                       results$bsample, method = 'PCA', k = 10)
visual <- DABB_visualize(Y, tot_read, gene_len, results$pweight,
                         results$bsample, method = 'ISOmap', k = 10)
write.table(visual, paste0('/Users/kexuanliang/Documents/singlecell/real_data/usoskin/visual.txt'), sep = '')

pc1 <- visual[,1]
pc2 <- visual[,2]

pc1.b1.t1  <- pc1 [ bat_ind == 1 & bio_ind == 1 ]
pc2.b1.t1  <- pc2 [ bat_ind == 1 & bio_ind == 1 ]
pc1.b1.t2  <- pc1 [ bat_ind == 1 & bio_ind == 2 ]
pc2.b1.t2  <- pc2 [ bat_ind == 1 & bio_ind == 2 ]
pc1.b1.t3  <- pc1 [ bat_ind == 1 & bio_ind == 3 ]
pc2.b1.t3  <- pc2 [ bat_ind == 1 & bio_ind == 3 ]
pc1.b1.t4  <- pc1 [ bat_ind == 1 & bio_ind == 4 ]
pc2.b1.t4  <- pc2 [ bat_ind == 1 & bio_ind == 4 ]

pc1.b2.t1  <- pc1 [ bat_ind == 2 & bio_ind == 1 ]
pc2.b2.t1  <- pc2 [ bat_ind == 2 & bio_ind == 1 ]
pc1.b2.t2  <- pc1 [ bat_ind == 2 & bio_ind == 2 ]
pc2.b2.t2  <- pc2 [ bat_ind == 2 & bio_ind == 2 ]
pc1.b2.t3  <- pc1 [ bat_ind == 2 & bio_ind == 3 ]
pc2.b2.t3  <- pc2 [ bat_ind == 2 & bio_ind == 3 ]
pc1.b2.t4  <- pc1 [ bat_ind == 2 & bio_ind == 4 ]
pc2.b2.t4  <- pc2 [ bat_ind == 2 & bio_ind == 4 ]

pc1.b3.t1  <- pc1 [ bat_ind == 3 & bio_ind == 1 ]
pc2.b3.t1  <- pc2 [ bat_ind == 3 & bio_ind == 1 ]
pc1.b3.t2  <- pc1 [ bat_ind == 3 & bio_ind == 2 ]
pc2.b3.t2  <- pc2 [ bat_ind == 3 & bio_ind == 2 ]
pc1.b3.t3  <- pc1 [ bat_ind == 3 & bio_ind == 3 ]
pc2.b3.t3  <- pc2 [ bat_ind == 3 & bio_ind == 3 ]
pc1.b3.t4  <- pc1 [ bat_ind == 3 & bio_ind == 4 ]
pc2.b3.t4  <- pc2 [ bat_ind == 3 & bio_ind == 4 ]

 p <- ggplot() + labs(title="Data Visualization", x="PC1", y="PC2")+
   theme(plot.title = element_text(hjust = 0.5))

 p <-  p + geom_point(data=data.frame(pc1.b1.t1 , pc2.b1.t1 ), aes(x=pc1.b1.t1 , y=pc2.b1.t1 , color="NF", shape='cold'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b1.t2 , pc2.b1.t2 ), aes(x=pc1.b1.t2 , y=pc2.b1.t2 , color="NP", shape='cold'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b1.t3 , pc2.b1.t3 ), aes(x=pc1.b1.t3 , y=pc2.b1.t3 , color="PHP", shape='cold'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b1.t4 , pc2.b1.t4 ), aes(x=pc1.b1.t4 , y=pc2.b1.t4 , color="TH", shape='cold'),    size=2)

 p <-  p + geom_point(data=data.frame(pc1.b2.t1 , pc2.b2.t1 ), aes(x=pc1.b2.t1 , y=pc2.b2.t1 , color="NF", shape='RT-1'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b2.t2 , pc2.b2.t2 ), aes(x=pc1.b2.t2 , y=pc2.b2.t2 , color="NP", shape='RT-1'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b2.t3 , pc2.b2.t3 ), aes(x=pc1.b2.t3 , y=pc2.b2.t3 , color="PHP", shape='RT-1'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b2.t4 , pc2.b2.t4 ), aes(x=pc1.b2.t4 , y=pc2.b2.t4 , color="TH", shape='RT-1'),    size=2)

 p <-  p + geom_point(data=data.frame(pc1.b3.t1 , pc2.b3.t1 ), aes(x=pc1.b3.t1 , y=pc2.b3.t1 , color="NF", shape='RT-2'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b3.t2 , pc2.b3.t2 ), aes(x=pc1.b3.t2 , y=pc2.b3.t2 , color="NP", shape='RT-2'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b3.t3 , pc2.b3.t3 ), aes(x=pc1.b3.t3 , y=pc2.b3.t3 , color="PHP", shape='RT-2'),    size=2)
 p <-  p + geom_point(data=data.frame(pc1.b3.t4 , pc2.b3.t4 ), aes(x=pc1.b3.t4 , y=pc2.b3.t4 , color="TH", shape='RT-2'),    size=2)


 ## visualize our result after quality control
 visual.qc <- DABB_visualize(Y[,-outliner$gene.num], tot_read[-outliner$gene.num], gene_len, results.qc$pweight,
                             results.qc$bsample, method = 'PCA', k = 10)
 visual.qc <- DABB_visualize(Y[,-outliner$gene.num], tot_read[-outliner$gene.num], gene_len, results.qc$pweight,
                             results.qc$bsample, method = 'ISOmap', k = 10)
pc1.qc <- visual.qc[,1]
pc2.qc <- visual.qc[,2]

# after qc
pc1.b1.t1.qc <- pc1.qc[ bat_ind == 1 & bio_ind == 1 ]
pc2.b1.t1.qc <- pc2.qc[ bat_ind == 1 & bio_ind == 1 ]
pc1.b1.t2.qc <- pc1.qc[ bat_ind == 1 & bio_ind == 2 ]
pc2.b1.t2.qc <- pc2.qc[ bat_ind == 1 & bio_ind == 2 ]
pc1.b1.t3.qc <- pc1.qc[ bat_ind == 1 & bio_ind == 3 ]
pc2.b1.t3.qc <- pc2.qc[ bat_ind == 1 & bio_ind == 3 ]
pc1.b1.t4.qc <- pc1.qc[ bat_ind == 1 & bio_ind == 4 ]
pc2.b1.t4.qc <- pc2.qc[ bat_ind == 1 & bio_ind == 4 ]

pc1.b2.t1.qc <- pc1.qc[ bat_ind == 2 & bio_ind == 1 ]
pc2.b2.t1.qc <- pc2.qc[ bat_ind == 2 & bio_ind == 1 ]
pc1.b2.t2.qc <- pc1.qc[ bat_ind == 2 & bio_ind == 2 ]
pc2.b2.t2.qc <- pc2.qc[ bat_ind == 2 & bio_ind == 2 ]
pc1.b2.t3.qc <- pc1.qc[ bat_ind == 2 & bio_ind == 3 ]
pc2.b2.t3.qc <- pc2.qc[ bat_ind == 2 & bio_ind == 3 ]
pc1.b2.t4.qc <- pc1.qc[ bat_ind == 2 & bio_ind == 4 ]
pc2.b2.t4.qc <- pc2.qc[ bat_ind == 2 & bio_ind == 4 ]

pc1.b3.t1.qc <- pc1.qc[ bat_ind == 3 & bio_ind == 1 ]
pc2.b3.t1.qc <- pc2.qc[ bat_ind == 3 & bio_ind == 1 ]
pc1.b3.t2.qc <- pc1.qc[ bat_ind == 3 & bio_ind == 2 ]
pc2.b3.t2.qc <- pc2.qc[ bat_ind == 3 & bio_ind == 2 ]
pc1.b3.t3.qc <- pc1.qc[ bat_ind == 3 & bio_ind == 3 ]
pc2.b3.t3.qc <- pc2.qc[ bat_ind == 3 & bio_ind == 3 ]
pc1.b3.t4.qc <- pc1.qc[ bat_ind == 3 & bio_ind == 4 ]
pc2.b3.t4.qc <- pc2.qc[ bat_ind == 3 & bio_ind == 4 ]

p1 <- ggplot()
 p1 <- p1 + geom_point(data=data.frame(pc1.b1.t1.qc, pc2.b1.t1.qc), aes(x=pc1.b1.t1.qc, y=pc2.b1.t1.qc, color="NF", shape='cold'),    size=2)
 p1 <- p1 + geom_point(data=data.frame(pc1.b1.t2.qc, pc2.b1.t2.qc), aes(x=pc1.b1.t2.qc, y=pc2.b1.t2.qc, color="NP", shape='cold'),    size=2)
 p1 <- p1 + geom_point(data=data.frame(pc1.b1.t3.qc, pc2.b1.t3.qc), aes(x=pc1.b1.t3.qc, y=pc2.b1.t3.qc, color="PHP", shape='cold'),    size=2)
 p1 <- p1 + geom_point(data=data.frame(pc1.b1.t4.qc, pc2.b1.t4.qc), aes(x=pc1.b1.t4.qc, y=pc2.b1.t4.qc, color="TH", shape='cold'),    size=2)

p1 <- p1 + geom_point(data=data.frame(pc1.b2.t1.qc, pc2.b2.t1.qc), aes(x=pc1.b2.t1.qc, y=pc2.b2.t1.qc, color="NF", shape='RT-1'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b2.t2.qc, pc2.b2.t2.qc), aes(x=pc1.b2.t2.qc, y=pc2.b2.t2.qc, color="NP", shape='RT-1'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b2.t3.qc, pc2.b2.t3.qc), aes(x=pc1.b2.t3.qc, y=pc2.b2.t3.qc, color="PHP", shape='RT-1'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b2.t4.qc, pc2.b2.t4.qc), aes(x=pc1.b2.t4.qc, y=pc2.b2.t4.qc, color="TH", shape='RT-1'),    size=2)

p1 <- p1 + geom_point(data=data.frame(pc1.b3.t1.qc, pc2.b3.t1.qc), aes(x=pc1.b3.t1.qc, y=pc2.b3.t1.qc, color="NF", shape='RT-2'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b3.t2.qc, pc2.b3.t2.qc), aes(x=pc1.b3.t2.qc, y=pc2.b3.t2.qc, color="NP", shape='RT-2'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b3.t3.qc, pc2.b3.t3.qc), aes(x=pc1.b3.t3.qc, y=pc2.b3.t3.qc, color="PHP", shape='RT-2'),    size=2)
p1 <- p1 + geom_point(data=data.frame(pc1.b3.t4.qc, pc2.b3.t4.qc), aes(x=pc1.b3.t4.qc, y=pc2.b3.t4.qc, color="TH", shape='RT-2'),    size=2)




####################################################################
#Visualization of the unfitted data.
#G_num <- length(Y[,1])
#C_num <- length(Y[1,])
#sample_num <- length(results$bsample[1,])
#visual <- DABB_visualize(Y, total.count, gene_length, matrix(1, G_num, C_num),
 #                        b_sample_mat = matrix(0, C_num, sample_num), method = 'PCA', k = 10)
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
