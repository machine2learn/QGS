### Install and load libs
if("genio" %in% rownames(installed.packages()) == FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("genio") # for read_plink()
}

if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
  install.packages("ggplot2") # for autoplot of prcomp
}

library("genio")
library("ggplot2")

### Read plink data, purely for the phenotypes
snps <- read_plink('hm3_hg38')

data <- snps$fam
data['Population'] <- data['pat']
data['y'] <- 1:1184

### Read QGS data and make plots

for (pop in c("ESN", "CHS", "FIN")) {

  # Load relevant QGS file
  qgs <- read.csv(file = sprintf('hm3_%s_ref.csv', pop), header = TRUE)
  data['dist'] <- colSums(qgs[, c(10:ncol(qgs))])


  # Barplot
  qgs_mean <- aggregate(data$dist, by=list(Category=data$Population), FUN=mean)
  qgs_std <- aggregate(data$dist, by=list(Category=data$Population), FUN=sd)
  bdata <- data.frame(Population = unique(data$Population), mean = qgs_mean$x, std = qgs_std$x)

  ggplot(bdata, aes(x=Population, y=mean)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2,
                   position=position_dodge(.9)) 

  ggsave(sprintf('plots/hm3_%s_bar.png', pop))

  # Scatterplot
  ggplot(data, aes(x=y, y=dist)) + 
    geom_point(aes(shape=Population, fill=Population)) + 
    scale_shape_manual(values=c(22,8,3,24,25,21,23,22,24,21,25)) + 
    scale_fill_manual(
      values=c(
        gray(0),gray(0),gray(0),gray(0.3),
        gray(0),gray(1), gray(0.8), gray(1),
        gray(1),gray(0.6),gray(1)
      )
    )

  ggsave(sprintf('plots/hm3_%s_scatter.png', pop))
}
