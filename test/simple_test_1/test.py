import pandas as pd
import numpy as np
import os
import sys
sys.path.append("..") # Find packages in parent dir

from qgs import qgs, write_vcf, fltcmp


# Sample matrices
S1 = np.array([[2, 1],[2, 1],[2, 1],[2, 1],[2, 1],[2, 1]]);

# Reference matrices
R1 = np.array([[0],[0],[0],[0],[0],[0]]);

# Calculate expected QGS
qgs_exp = qgs(S1, R1)

# Write vcf files
write_vcf("sample.vcf", S1)
write_vcf("reference.vcf", R1)

# Run external QGS
os.system("../../qgs --sample sample.vcf --reference reference.vcf --genes ../testgenes.gtf --out test.csv --maf 0")

# Read QGS output file
qgs_values = pd.read_csv("test.csv", header=0)
qgs_gene1 = qgs_values.loc[0][9:].to_numpy().astype(float)

# Test
if np.sum(np.abs(np.subtract(qgs_exp, qgs_gene1))) < 0.000001:
  print("Test passed.")
  sys.exit(0)

print("Test FAILED.")
sys.exit(1)