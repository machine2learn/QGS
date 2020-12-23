import pandas as pd
import numpy as np
import os
import sys

sys.path.append("..")  # Find packages in parent dir

from qgs import qgs, write_vcf, fltcmp

# number_samples = 100
# number_variants = 15
failed_tests = 0
passed_tests = 0

for number_samples, number_variants in zip([10, 10, 100, 100, 100, 1000], [2, 2, 10, 10, 15, 15]):
    # Sample matrices
    S1 = np.random.uniform(0, 2, number_variants * number_samples).reshape(number_variants, number_samples)

    # Reference matrices
    R1 = np.random.uniform(0, 2, number_variants * number_samples).reshape(number_variants, number_samples)

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
    if np.sum(np.abs(np.subtract(qgs_exp, qgs_gene1))) < 0.0003:
        print("Test passed.")
        passed_tests += 1
    else:
        print(np.sum(np.abs(np.subtract(qgs_exp, qgs_gene1))))
        print("Test FAILED.")
        failed_tests += 1

print(f"-----------------------------------------------------------")
print(f"{passed_tests} tests passed and {failed_tests} tests failed")

if failed_tests >= 1:
    sys.exit(1)
else:
    sys.exit(0)
