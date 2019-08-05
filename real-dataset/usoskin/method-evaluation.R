## PCA of raw data
raw <- exprs(eset)

totalcount = function (ei)
{
        sums = colSums(ei)
        eo = t(t(ei)*mean(sums)/sums)
        return(eo)
}

tc <- totalcount(raw)

library(rARPACK)
fastpca <- function(expr, scale=FALSE) {
        svd_raw <- svds(scale(t(expr), center=TRUE, scale=scale), k=3, nu=3, nv=0)
        pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:3])
        return(pc_raw)
}

vars <- rowVars(log1p(tc))
names(vars) <- rownames(tc)
vars <- sort(vars, decreasing = TRUE)
vargenes <- names(vars)[1:1000]

pc_raw <- fastpca(log1p(raw[vargenes,]))
pc_tc <- fastpca(log1p(tc[vargenes,]))

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

level1 <- as.factor(pData(eset)$`Picking sessions`)
level2 <- as.factor(pData(eset)$`Level 1`)

colMerged <- col1[level1]
colCl <- col2[level2]

par(mar=c(4,4,2,6))
par(mgp=c(2,0.5,0))
plot(pc_raw, col=colMerged, pch=20, main="PCA RAW", cex=.5)
usr <- par("usr")
x <- usr[2]*1.05
y <- usr[4]*0.2
legend(x,y, levels(level1), fill=col1, xpd=T)

par(mar=c(4,4,2,6))
par(mgp=c(2,0.5,0))
plot(pc_raw, col=colCl, pch=20, main="PCA RAW", cex=.5)
usr <- par("usr")
x <- usr[2]*1.05
y <- usr[4]*0.2
legend(x,y,levels(level2), fill=col2, cex=0.5, xpd=T)
