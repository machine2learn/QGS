import numpy as np

def fltcmp(lhs, rhs):
  """Simple helper function to round and compare two floats

  Args:
    lhs (float)
    rhs (float)
  
  Both arguments expect subjects in rows and variants in columns

  Returns:
    int: 0 if equal; -1 if lhs < rhs; +1 otherwise

  """
  
  if round(lhs, 5) == round(rhs, 5):
    return 0
  elif lhs < rhs:
    return -1
  return 1

def qgs_single(s, R):
  """Calculate QGS score

  Args:
    s (numpy.array): Vector containing single sample dosages
    R (numpy.array): Vector or matrix containing reference dosages
  
  Both arguments expect subjects in rows and variants in columns

  Returns:
    float: The QGS sum score of the inputs for one individual

  """

  return np.sum(np.abs(np.subtract(s, R))) / (2 * np.product(np.shape(R)))


def qgs(S, R):
  """Calculate QGS score

  Args:
    S (numpy.array): Vector or matrix containing sample dosages
    R (numpy.array): Vector or matrix containing reference dosages
  
  Both arguments expect subjects in rows and variants in columns

  Returns:
    np.array of float: The QGS sum score for all samples

  """
  
  if np.ndim(S) == 1:
    return qgs_single(S, R)
  
  return np.apply_along_axis(qgs_single, axis=1, arr=S, R=R)
  
def write_vcf(filename, D, startpos = 1):
  """Write simple VCF file
  
  Function has a side effect: it keeps its own internal counter

  Args:
    filename (string): Path to output file
    D (numpy.array): Vector or matrix containing dosages
    startpos (unsigned int): Output base pair start position
  
  D expects subjects in rows and variants in columns

  Returns:
    boolean: True on success or False otherwise

  """
  
  out = """##fileformat=VCFv4.1
##For format description see https://samtools.github.io/hts-specs/VCFv4.1.pdf
##FORMAT=<ID=DS,Number=1,Type=String,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

  f = open(filename, "w")
  f.write(out)
  
  n = len(D) # number of samples
  if D.ndim == 1:
    n = 1
  
  d = len(D) # number of variants
  if D.ndim > 1:
    d = D.shape[1]
  
  for i in range(n):
    f.write("	SubID%s" % (i))
  
  f.write("\n")
  
  for idx, col in enumerate(D.T):
    f.write("1\t%d\t1:%d\tA\tG\t100\tPASS\tVT=SNP\tDS\t" % (startpos + idx, startpos + idx))
    if n == 1:
      f.write(str(col))
    else:
      f.write('\t'.join([str(ds) for ds in col]))
    f.write("\n")

  f.close()


