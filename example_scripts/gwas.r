# Read the file name from the command line argument
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Provide input file name")
}
qgs_file = args[1]

# Read the QGS file
qgs_values <- read.csv(file = qgs_file, header=T)
rows <- nrow(qgs_values)
cols = ncol(qgs_values)

# Remember that the first 9 values from any row are meta-information
# (gene name, location, etc) and the actual QGS values start at idx 10
nsubjects <- cols - 9

# Generate a deterministic phenotype for the gwas
pheno <- seq(0, 1, length.out = nsubjects)

# Loop over every row/gene
for (gene_idx in 1:rows) {

  # Extract the QGS values for this gene
  qgs_row <- t(as.matrix(qgs_values[gene_idx, c(10:cols)]))
  
  # Perform linear regression
  mdl = lm(pheno ~ qgs_row)
  
  # Print output
  cat(sprintf("%s\t%.2f\t%e\n",
    qgs_values[gene_idx, 1], # gene name
    summary(mdl)$coefficients[2,1], # coefficient/beta
    summary(mdl)$coefficients[2,4] # p-value
  ))

}