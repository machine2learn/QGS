# Principal component analysis comparison

This example compares SNP-based PCA to QGS-based PCA. It uses public HapMap3 data to demonstrate that the first two principal components from SNPs and QGS are close to identical (a correlation of r>0.99). In addition, it creates two population stratification plots based on the two PCs. This analysis is used in the QGS publication.

The Linux shell script `hm3.sh` downloads and extracts the HapMap3 data, converts the data to Plink binary format using the `plink` command line tool, calculates QGS values, and calls R to do the PCA.

The R script `hm3.R` calculates PCs for both the SNPs and QGS and plots these.

The two pdf files contain the output of the script: two population stratification plots.
