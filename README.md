# QGS: Quantitative Genetic Scoring

QGS is a computer program to create genetic variables from SNP data. For every individual, it calculates the genetic distance to a reference population. QGS works with both VCF or PLINK format data. The pre-print of the QGS manuscript is available at https://doi.org/10.1101/2020.12.15.422886.

# Quickstart

Executable files for generic Linux and Windows x64 are available. The v1.0 generic Linux version can be [downloaded here](https://github.com/machine2learn/QGS/releases/download/v1.0/qgs-x86_64), and the Windows version can be [downloaded here](https://github.com/machine2learn/QGS/releases/download/v1.0/qgs-win64.exe). More recent binaries may be available from the repository. **This document assumes that the QGS executable is named "qgs" (or "qgs.exe" for Windows).** If you use the precompiled executables, it is recommended to rename them from "qgs-x86_64" to "qgs" (or "qgs-win64.exe" to "qgs.exe" for Windows).

If you want to compile the program yourself, see **Compiling** below.

Run the program from the command line using `./qgs`. The command `./qgs --help` will display all arguments.

There are four mandatory arguments, each of which is the path to a file:

    ./qgs --sample samplefile.vcf --reference referencefile.vcf --genes genefile.gtk --out outputfile.csv
    
Sample and reference files can be in VCF or PLINK bed format and will be treated as such based on file extension. All input and output files may be gzip compressed, which is also auto-detected based on extension.

This command will generate an `outputfile.csv` file containing the QGS values.

# Compiling

To compile the executable, you'll need `make`, libz, and a C++11-compatible compiler. On most Linux systems, these should be available. The Makefile assumes `g++` as the compiler. If you use another compiler, edit Makefile line 4.

    git clone https://github.com/machine2learn/QGS.git
    cd QGS
    make
    ./qgs --version

# Command line arguments

### Required arguments:

Option|Description|Example
------|-----|------
--sample [path]|Path to sample input file in VCF or Plink format. May be gz-compressed. (1) VCF files can be processed in either genotype (GT) or dosage (DS) format. (2) Plink binary bed/bim/fam files can also be processed. All three files must be present in the same directory. The path given to --sample must be that of a single .bed file, e.g. `--sample plinkfile.bed`. The other two files (bim/fam) will be found automatically. (3) Plink dosage files (such as produced by the [RICOPILI](https://doi.org/10.1093/bioinformatics/btz633) pipeline) are also supported. The path(s) given to --sample must be one or more `.dosage(.gz)` files. The reading from multiple Plink dosage files is supported, including with a wildcard. For example `--sample chr22_*.dosage.gz` will find `chr22_000_020.dosage.gz`, `chr22_020_040.dosage.gz` etc. Similarly, `--sample *.dosage.gz` will use all `dosage.gz` files in the directory. Corresponding `.map` files must be available in the same directory. | `--sample samplefile.vcf.gz`
--reference [path]|Path to reference input file in VCF or Plink format. See `--sample` for details. For best results, make sure the reference panel matches the sample population by geographical region. If the sample file is imputed, the imputation reference panel is a good choice. The sample file itself can also be used as a reference.| `--reference referencefile.vcf.gz`
--genes [path]|Path to gene location annotation in Gene Transfer Format. my be gz-compressed. We recommend using the Gencode releases.| `--genes gencode.basic.annotation.gtf.gz`
--out [path]|filename for writable output file. any existing file will be overwritten without warning.|`--out outputfile.csv`

**Make sure that the sample file, the reference file, and the gene locations share the same mapping/build** (e.g. Grch38/hg38 or b37/hg19). Mismatches between the sample and reference will likely be detected, but mismatches between the genotypes and the gene annotation cannot be detected.

### Optional arguments

###### Misc

Flag|Description|Example
------|-----|------
--help|Provides information about command line arguments|`--help`
--version|Prints QGS version and exits. Other arguments are ignored.|`--version`


###### Filtering

Flag|Description|Example
------|-----|------
--chr [value]|Only calculate QGS for listed chromosome. Value must be numeric. X = 23, Y = 24, MT = 25.|`--chr 22`
--maf [value]|Minor allele frequency threshold. Any variants with a lower observed MAF in either sample or reference will be excluded. Default value: 0.01|--maf 0.05
--gtf-filter [key]=[value] [[key2=value2]]|Key-value pairs for only including specific genes/transcripts or types. Multiple key-value pairs can be provided. Including the option `type=gene` is recommended to exclude transcripts, exons, etc. The `type=` argument filters specifically based on column #3 in the GTF file. Other filter options are based on column #9 from the GTF file and can take any key/value combination, e.g. `gene_type=protein_coding` (excludes non-protein-coding genes) or `gene_name=BRK1` (only calculate QGS for BRK1 gene).|`--gtf-filter type=gene gene_type=protein_coding`
--include-snps [path]|Path to a file containing variant names to include. Other variants will be excluded. Variant-names can be rs-numbers (if those are provided in the genetic input files) or chr:location. Variants must be separated by white space, e.g. one on every line. Cannot be combined with `--exclude-snps`| `--include-snps list_of_snps.txt`
--exclude-snps [path]|Path to a file containing variant names to exclude. Variant-names can be rs-numbers (if those are provided in the genetic input files) or chr:location. Variants must be separated by white space, e.g. one on every line. Cannot be combined with `--include-snps`| `--exclude-snps list_of_snps.txt`

###### Command line verbosity

Flag|Description|Example
------|-----|------
--verbose|Increase the information that QGS prints to the screen. Extra information includes number of individuals read and which scores were outputted.|`--verbose`
--debug|Increase the information that QGS prints to the screen even more. Extra information includes every successful read of a variant and which variants are includes and excluded for every score .|`--debug`
--trace|Print every decision that QGS makes. This will print **a lot** of output to the screen and significantly slow down the program. Only use this option if you need to find out *why* QGS is making some decision and `--debug` is not telling you enough.|`--trace`

If multiple verbosity options are provided, QGS will use the most verbose flag provided.

###### Flanking regions

By default, no flanking regions are included.

Flag|Description|Example
------|-----|------
--flank [value]|Symmetrical flanking region to be included for every QGS region. Value should be a positive whole number representing the flanking region in kb. Default value: 0.|`--flank 5`
--pre-flank [value]|Asymmetrical flanking region. Specifies the number of kb to be included in front of every QGS region. Regions with reverse strand orientation are handled automatically, provided the information is available in the GTK file. If both `--flank` and `--pre-flank` are provided, `--pre-flank` takes priority.|`--pre-flank 15`
--post-flank [value]|Asymmetrical flanking region. Specifies the number of kb to be included after every QGS region. Regions with reverse strand orientation are handled automatically, provided the information is available in the GTK file. If both `--flank` and `--post-flank` are provided, `--post-flank` takes priority.|`--post-flank 35`

###### Missing data

By default, any variant that has missing data for one or more subjects is excluded. In other words, a variant must be available for all subjects to be included in QGS. The flags below can change this behaviour **for the sample file**. Variants with missing information in the reference population will always be excluded.

Flag|Description|Example
------|-----|------
--fill-missings|Crudely impute any missing data points, assuming a 0/0 genotype.|`--fill-missings`
--allow-missings|Do not exclude variants with missing data, but instead calculate QGS values for everyone with complete data. For subjects with incomplete data, print 'NaN' (not a number) for that particular value.|`--allow-missings`

###### Input data handling

Flag|Description|Example
------|-----|------
--hard-calls|Force QGS to use hard calls. By default, QGS automatically chooses to use dosage information if available. With this flag, QGS always uses hard calls.|`--hard-calls`
--weight-by [label]|Provide a weight for each variant. Only works with VCF input file (sample or reference). By default, QGS values are calculated without variant weights, i.e. every variant is counted equally. With this option, it is possible to provide a VCF INFO field label which contains the weight for each variant. If provided, variants without a weight from the VCF file are ignored.|`--weight-by R2`

###### Output file options

Flag|Description|Example
------|-----|------
--delimiter [char]|By default, the output file is a comma (',') separated values file. Any 1 character provided here will replace the comma as delimiter in the output file.|`--delimiter ";"`
--output-variants|By default, the QGS output file contains the number of variants included used for each value. With this flag, that number is replaced by a pipe ('\|') delimited list of variant names. This helps make visible which specific variants are included in any QGS value.|`--output-variants`

# Output file format

The output file produced by QGS is a delimited comma-separated values (csv) file. It can be opened with a spreadsheet application (Excel, Calc, etc) or text editor. The format is as follows:

gene_name|gene_id|chr|start|stop|Nsample|Nref|num_loci|total_num_loci|SubjectID1|SubjectID2|...
---------|-------|---|-----|----|-------|----|--------|--------------|----------|----------|------
BRK1|ENSG00000254999.4_5|3|10157359|10168874|5|2504|55|470|0.178282|0.21087|...
...|||||||||||...

The output file has (9 + the number of samples) columns and each row contains one genetic location with QGS values for every sample individual. The columns contain the following information:

1. **gene_name**: Gene name (string) as given in the GTF file by the "gene_name" attribute. If none is provided, it will be auto-generated as chr:start-stop (e.g. "22:1000-2000").
2. **gene_id**: Gene ID (string) as given in the GTF file by the "gene_id" attribute. If none is provided, it will be identical to the gene_name.
3. **chr**: Chromosome (integer) of the genetic region.
4. **start**: Start location (integer) in base pairs of the genetic region, including any flanking region.
5. **start**: Stop location (integer) in base pairs of the genetic region, including any flanking region.
6. **Nsample**: Number (integer) of individuals for which the QGS values are available.
7. **Nref**: Number (integer) of individuals that were present in the reference population.
8. **num_loci**: Either (a) Number (integer) of variants included in the current value or (b) a pipe-separated list (string) of variant names included in the current value. Contents in determined by the `--output-variants` flag.
9. **total_num_loci**: Number (integer) of available variants that were observed to lie inside the current genetic region.
10. And onwards: QGS value for the current subject (column) for the current genetic region.

Example scripts to run a GWAS in both Python and R are available in the `example/` directory.

# Gene annotation files

The requires `--genes` command line option needs a gene annotation file in GTF -- Gene Transfer Format. A description of this format can be found here: https://www.ensembl.org/info/website/upload/gff.html

For our paper focussing on human genetics, we have used the annotations provided by [GENCODE](https://www.gencodegenes.org/). These can be downloaded from here: [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)

Make sure to select a gene annotation that matches the [build/reference assembly](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19) of your genetic data.

# Tutorial

In this tutorial, we will use the 1000 Genomes data to calculate QGS values for all protein coding genes on chromosome 22 and run a gene-based GWAS.

This tutorial assumes that you have a Linux command line (bash) environment.

### Setup

First, we will download and compile the QGS program:

    git clone https://github.com/machine2learn/QGS.git
    cd QGS
    make
    ./qgs --version
    
This should print the current QGS version. If you run into problems here, see the **Compile** section above.

Next, we will download the genetic data and gene annotation files (around 250MB in total):

    mkdir data
    cd data
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.basic.annotation.gtf.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz

Always make sure that the genetic data files and the gene annotation file have the same mapping/build, in this case GRCh38.

### Running QGS

Now we will run QGS and create an output file "out.csv" containing QGS gene values:

    cd ..
    ./qgs --sample data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --reference data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --genes data/gencode.v35.basic.annotation.gtf.gz --out out.csv --gtf-filter type=gene gene_type=protein_coding
    
This command will take a few minutes to run. You may see a few warnings about skipped variants.

In this example, we use the 1000 Genomes data as both the `--sample` file and the `--reference` file. The location of the gene annotation file is given with `--genes`. Finally, by default QGS constructs values for all entries in the annotation file. However, that file also contains locations of transcripts, exons, and other things that are not protein coding genes. Therefore, we include the filter `--gtf-filter type=gene gene_type=protein_coding` to only select protein coding genes.

After the run, you will see a short summary of how many variants were read from the input file and how many were actually used:

    [INFO] Sample statistics:
      Loci read: 1100291
      Loci used: 95918 (8.72%)
    Run completed

In the run above, about 9% of all variants from the input file were actually used in QGS values. This may seem low, but is not unexpected: remember that genes (and especially protein coding genes) do not cover our entire genome. In fact, most genetic information lies outside of gene regions. In addition, a lot of the variants in the 1000 Genomes data are below the default MAF threshold of 0.01.

This summary information is important, because you can quickly check if something went wrong. For example, if you have a mismatch between sample and reference mapping, the "loci used" will be very low.

With the above command, we have created the "out.csv" QGS file. This file contains the generated genetic variables. The format and contents of the file is detailled above in the section **Output file format**. These variables can be used in your followup analyses.

### Example of GWAS analysis

To show how the generated QGS may be used, example code to perform a GWAS analysis in R and Python has been provided in the `example_scripts/` directory.

    cd example_scripts/
    

Then a GWAS analysis in Python can be started with the command 

    python3 gwas.py ../out.csv
    
If you are missing Python libraries, make sure to install `pandas`, `numpy`, and `scipy`. In most environments, this can be done with the `pip3 install pandas numpy scipy` command.

Similarly, the R GWAS script can be started like this:

    Rscript gwas.r ../out.csv

You only need to run one these commands, because both do the same thing but in different programming languages.

The script must be provided with the location of the "out.csv" file as seen above. The output of the script is a list consisting of a gene name, a regression coefficient/beta, and a p-value for that gene. If you're using the "out.csv" file generated by the tutorial commands above, the first three lines of output should be something like this:

    OR11H1	-0.12	7.203238e-02
    POTEH	0.64	1.062172e-02
    CCT8L2	0.05	4.425209e-01
    
We can check what the top 5 most significantly associated genes are by sorting by the third column (assuming a Linux Bash-like command line):

    python3 gwas.py ../out.csv | sort -gk3 | head -5
    Rscript gwas.r ../out.csv | sort -gk3 | head -5
    
The output should be something like

    BRD1	1.64	3.047018e-22
    GNB1L	2.13	2.448744e-20
    SNRPD3	0.93	4.791433e-20
    ZBED4	0.94	8.238983e-20
    ALG12	1.01	6.859289e-19

Keep in mind that these results are meaningless, because the example scripts generate a fake deterministic phenotype for the GWAS regression. For your own analyses, you will need to edit the scripts to include your own phenotype (and likely covariates), or run the analyses using any other software you prefer (e.g. SPSS, MATLAB, etc).

Unless you have a reason not to, we recommend including QGS-based principal components as covariates in your association analyses. For the analyses in the QGS paper, we included 10 to 20 PCs.

### PCA analysis

See here for the PCA analysis as used in the QGS paper: [example_scripts/pca_comparison/](example_scripts/pca_comparison/).

### Different command line arguments

By changing some default settings, we can choose to include more variants in the QGS scores.

For example, we can lower the MAF threshold from the default 0.01 to 0, effectively removing the MAF threshold by including the `--maf 0` flag. Note that we also changed the output file with `--out out_maf0.csv`.

    ./qgs --sample data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --reference data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --genes data/gencode.v35.basic.annotation.gtf.gz --out out_maf0.csv --gtf-filter type=gene gene_type=protein_coding --maf 0
    
This results in the following summary statistics:

    [INFO] Sample statistics:
      Loci read: 1100291
      Loci used: 600364 (54.6%)
    Run completed

From this, we can see that many of the variants in the 1000 Genomes data are rare variants (<1% MAF).

We can also include symmetrical (with `--flank`) or asymmetrical flanking regions:

    ./qgs --sample data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --reference data/ALL.chr22_GRCh38.genotypes.20170504.vcf.gz --genes data/gencode.v35.basic.annotation.gtf.gz --out out_flank.csv --gtf-filter type=gene gene_type=protein_coding --pre-flank 15 --post-flank 35

This results in

    [INFO] Sample statistics:
      Loci read: 1100291
      Loci used: 136147 (12.4%)
    Run completed

For a full overview of possible arguments, see the section **Command line arguments** above.

# Funding

This project has received funding from the European Union’s Seventh Framework Programme for research, technological development and demonstration under grant agreement no 602805 - AGGRESSOTYPE. This work reflects only the author’s views and the European Union is not liable for any use that may be made of the information contained therein. This work is part of the research programme Computing Time National Computing Facilities Processing Round pilots 2018 with project number 17666, which is (partly) financed by the Dutch Research Council (NWO). This work was partially carried out on the Dutch national e-infrastructure with the support of SURF Cooperative.
