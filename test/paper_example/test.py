import pandas as pd
import numpy as np
import os
import sys
from ..qgs import qgs, write_vcf, fltcmp


# Sample matrix
S = np.array([0.2, 1.2]);

# Reference matrix
R = np.array([[0.5, 0.8], [1.1, 0.1], [1.4, 0.3]])

# Write vcf files
write_vcf("sample.vcf", S)
write_vcf("reference.vcf", R)

# Run external QGS
os.system("../../qgs --sample sample.vcf --reference reference.vcf --genes ../testgenes.gtf --out test.csv --maf 0")

# Read QGS output file
qgs_values = pd.read_csv("test.csv", header=0)
qgs_value = qgs_values.SubID0[0]

# Test
if fltcmp(qgs(S, R), qgs_value) == 0:
  print("Test passed.")
  sys.exit(0)

print("Test FAILED.")
sys.exit(1)