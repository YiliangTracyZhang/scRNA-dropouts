##################
##pre-processing##
##################
# as the read count data does not provide gene length information, so we provide three ways to map gene id to gene information
# route 1:
# library (EDASeq)
# gene_len_list <- getGeneLengthAndGCContent(ensembl_list, "hsa")
# gene_len_list <- getENSG(ensembl_list, mart.gene)
# route 2
library(biomaRt)
mart.gene = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
getENSG <- function(gene, mart = mart.gene) {
results <- getBM(attributes = c("ensembl_gene_id","start_position","end_position","transcript_length","hgnc_symbol"),
                filters    = "hgnc_symbol", values = gene, mart = mart)
return(results$transcript_length)
}

# route 3:
# getting human gene length data in advance. as the former methods don't perform so well, we adopt to this method.
# linux: wget ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
#
blakeley <- read.table('/Users/kexuanliang/Documents/singlecell/git/scRNA-dropouts/real-dataset/blakeley.txt', sep = '')
rownames(blakeley) <- blakeley[,1]
cell_info <- blakeley[1,-1]
gene_info <- blakeley[-1,1]
C_num <- length(cell_info)
G_num <- length(gene_info)
expr.bkl <- as.matrix(blakeley[-1,-1])
expr.bkl <- matrix(as.numeric(expr.bkl), G_num, C_num)

ensembl_list <- as.character(gene_info)
ref <- read.table('/Users/kexuanliang/Documents/singlecell/real_data/hg19/genes.tsv', sep ='', header = F)
names(ref) <- (ensembl, gene.id)
gene_len <- ref[ match(ensembl_list, ref$ensemble), 2]
gene_len <- getENSG(gene_len, mart.gene)
#gene_len_list <- getENSG(ensembl_list, mart.gene)
# library (EDASeq)
# gene_len_list <- getGeneLengthAndGCContent(ensembl_list, "hsa")

