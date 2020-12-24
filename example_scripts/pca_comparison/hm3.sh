# Download & extract HM3 data
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2
bunzip2 hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2
bunzip2 hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2

# Convert to Plink binary format exclusing non-autosomal and missing
# This command below uses Plink v1.90
plink --file hapmap3_r2_b36_fwd.consensus.qc.poly --make-bed --autosome --geno 0 --out hm3

# Include HM3 population names in Plink fam file
# Because Plink phenotypes can only be numeric, we abuse the 
# paternal ID column (#3) of the fam file to store the population names
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_2/relationships_w_pops_121708.txt
mv hm3.fam hm3.fam.org
awk 'FNR==NR { a[$1]=$7; next; } { $3 = a[$1]; print $0; }' relationships_w_pops_121708.txt hm3.fam.org > hm3.fam

# Run QGS
# HM3 build is NCBI36, we download a matching gene annotation file:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_3c/gencode.v3c.annotation.NCBI36.gtf.gz
qgs --sample hm3.bed --reference hm3.bed --genes gencode.v3c.annotation.NCBI36.gtf.gz --out qgs_hm3.csv --maf 0 --gtf-filter type=gene

### Run PCA and make plots in R
Rscript hm3.R