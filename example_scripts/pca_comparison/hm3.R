### Install and load libs
if("genio" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("genio") # for read_plink()
}

if("ggfortify" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggfortify") # for autoplot of prcomp
}

library("genio")
library("ggfortify")

### SNP PCA
### This takes a lot of memory (>20GB)
### Better to use alternatives such as Plink, GCTA64, flashpcaR, etc
### We use pccomp here for 100% compatibility with QGS PCA
snps <- read_plink('hm3')
spca <- prcomp(t(snps$X), scale = TRUE, rank = 2, retx = TRUE)

### QGS PCA
qgs <- read.csv(file = 'qgs_hm3.csv', header=T)
qpca <- prcomp(t(qgs[, c(10:ncol(qgs))]), scale = TRUE, rank = 2, retx = TRUE)

### The PCs have a high correlation:
message(sprintf("Correlation PC1: %.4f\n", cor(spca$x[,'PC1'], qpca$x[,'PC1'])))
message(sprintf("Correlation PC2: %.4f\n", cor(spca$x[,'PC2'], qpca$x[,'PC2'])))

### PCA gives arbitrary sign, so transform data to match each other for plot
if (cor(spca$x[,'PC1'], qpca$x[,'PC1']) < 0)
  spca$x[,'PC1'] = spca$x[,'PC1'] * -1
if (cor(spca$x[,'PC2'], qpca$x[,'PC2']) < 0)
  spca$x[,'PC2'] = spca$x[,'PC2'] * -1
  
### Grab the population labels from the Plink paternal id column
### and rename the column for correct plot legend name
poplabels = snps$fam
poplabels['Population'] = poplabels['pat']

### Make QGS plot
autoplot(qpca, data = poplabels, shape = 'Population', fill = 'Population') +
  scale_shape_manual(values=c(22,8,3,24,25,21,23,22,24,21,25)) +
  scale_fill_manual(
    values=c(
      gray(0),gray(0),gray(0),gray(0.3),
      gray(0),gray(1), gray(0.8), gray(1),
      gray(1),gray(0.6),gray(1)
    )
  )
ggsave("pca_qgs.pdf")

### Make SNP plot
autoplot(spca, data = poplabels, shape = 'Population', fill = 'Population') +
  scale_shape_manual(values=c(22,8,3,24,25,21,23,22,24,21,25)) +
  scale_fill_manual(
    values=c(
      gray(0),gray(0),gray(0),gray(0.3),
      gray(0),gray(1), gray(0.8), gray(1),
      gray(1),gray(0.6),gray(1)
    )
  )
ggsave("pca_snp.pdf")
