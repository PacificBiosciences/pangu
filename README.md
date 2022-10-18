# Pangu: A star-typer for long-read PGx applications

The PacBio CYP2D6 star-typer (Pangu) for whole genome sequencing, capture, enrichment, and amplicon data. This caller leverages the accuracy and length of HiFi reads to confidently call star alleles, including structurally diverse alleles. The caller produces [pharmvar](https://www.pharmvar.org/gene/CYP2D6) compatible definitions. 

# Demo Dataset

See `demo_data` directory for example HiFi datasets.

## Installation

### Dependencies
 - [pysam](https://github.com/pysam-developers/pysam)
 - [mappy](https://pypi.org/project/mappy/)
 - [pandas](https://pandas.pydata.org/)
 - [numpy](https://numpy.org/)
 - [scipy](https://scipy.org/)
 - [scikit-learn](https://scikit-learn.org/stable/index.html)

### Installing in a Conda enviroment 
```
git clone https://github.com/PacificBiosciences/pangu.git
conda create --name cyp python=3.9
cd pangu
python setup.py install
```

### Reference

You need to supply your own GRCh38 reference (with .fai index), as this is too large for the repo.  Use just chromosome 22 to keep analysis faster.

```
samtools faidx GRCh38_full.fasta chr22 > GRCh38_chr22.fasta && samtools faidx GRCh38_chr22.fasta
```

Then add the full path to the reference in the file `genes/CYP2D6.yaml` (first setting at the top).

## To run the caller

```
# input bam contains HiFi reads aligned to GRCh38
./pangu -p outdir/sample_name --verbose <inBam>

# Restrict variants used for phasing of reads to positions contained in a vcf (e.g. from DeepVariant).
./pangu -p outdir/sample_name --verbose --vcf <inVcf> <inBam>

# Run on multiple samples 
./pangu -p outdir/sample_name --verbose path/to/data/*bam
# or
./pangu -p outdir/sample_name --verbose --bamFofn bams.fofn

# Export a bam containing just the reads analyzed for CYP2D6.  
# Reads are labeled by allele and colored by haplotype.
# (use the --grayscale option if you don't like the colors!)
# use logLevel DEBUG for more details on the allele calling
./pangu -p outdir/sample_name --verbose -x --logLevel DEBUG  <inBam>  
```

## Caller Outputs

Results will be generated with the prefix given, else in the cwd.  Diplotype calls are found in the \<prefix\>[/\_]call.log file.
If you include the `-v,--verbose` option, results will be printed to the screen.
```
./pangu -p demo_data/example/NA20129 -x -v demo_data/NA20129_wgs.GRCh38.bam
2022-10-18 15:27:27,861 | INFO | Processing demo_data/NA20129_wgs.GRCh38.bam
2022-10-18 15:27:36,423 | INFO | Identifying SV signatures
2022-10-18 15:27:39,752 | INFO | Phasing and Calling Alleles
2022-10-18 15:27:52,644 | INFO | Calling Diplotype
2022-10-18 15:27:52,645 | INFO | Diplotype called: CYP2D6 *4x2/*5
2022-10-18 15:27:52,645 | INFO | Exporting demo_data/example/NA20129.labeledReads.bam
```
Color:

![demo image](images/NA20129_igv_color.png?raw=true "Labeled Reads")

Grayscale:

![demo image](images/NA20129_igv_grayscale.png?raw=true "Labeled Reads (grayscale)")

## Release History
* 0.1.0 - Initial Release (10/18/2022)


## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

FOR RESEARCH USE ONLY. NOT FOR USE IN DIAGNOSTICS PROCEDURES.
