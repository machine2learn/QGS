import sys
import pandas as pd
import numpy as np
from scipy import stats

# Read the file name from the command line argument
if len(sys.argv) != 2:
  sys.exit("Provide input file name")
qgs_file = sys.argv[1]

# Read the QGS file
qgs_values = pd.read_csv(qgs_file, header=0)

# Remember that the first 9 values from any row are meta-information
# (gene name, location, etc) and the actual QGS values start at idx 10
nsubjects = qgs_values.shape[1] - 9

# Generate a deterministic phenotype for the gwas
pheno = np.linspace(0, 1, nsubjects)

# Loop over every row/gene
for gene_idx in qgs_values.index:
  # Extract the QGS values for this gene
  qgs_row = qgs_values.loc[gene_idx][9:].to_numpy().astype(float)

  # Perform linear regression
  beta, intercept, r_value, p_value, std_err = stats.linregress(qgs_row, pheno)
  
  # Print output
  print("%s\t%.2f\t%e" % (
    qgs_values.loc[gene_idx][0], # gene name
    beta, # coefficient/beta
    p_value # p value
  ))

