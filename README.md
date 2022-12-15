# Pangu: A star-typer for long-read PGx applications

The PacBio star-typer (Pangu) for whole genome sequencing, capture, enrichment, and amplicon data. This caller leverages the accuracy and length of HiFi reads to confidently call star alleles, including structurally diverse alleles. The caller produces [pharmvar](https://www.pharmvar.org/gene/CYP2D6) compatible definitions. Currently calling CYP2D6.

## Demo Dataset

See `demo_data` directory for example HiFi datasets.

## Installation

### Dependencies
 - python >=3.9
 - [pysam](https://github.com/pysam-developers/pysam)
 - [mappy](https://pypi.org/project/mappy/)
 - [pandas](https://pandas.pydata.org/)
 - [numpy](https://numpy.org/)
 - [scipy](https://scipy.org/)
 - [scikit-learn](https://scikit-learn.org/stable/index.html)

### Reference
GRCh38 chr22 is included in the repo for convenience.  

### From github 
```
git clone https://github.com/PacificBiosciences/pangu.git
cd pangu
pip install .
```

### From bioconda (please wait for bioconda PR to be merged)
```
conda install -c bioconda pangu
```

## To run the caller

```
# input bam contains HiFi reads aligned to GRCh38
pangu -p outdir/sample_name --verbose <inBam>

# Restrict variants used for phasing of reads to positions contained in a vcf (e.g. from DeepVariant).
pangu -p outdir/sample_name --verbose --vcf <inVcf> <inBam>

# Run on multiple samples 
pangu -p outdir/sample_name --verbose path/to/data/*bam
# or
pangu -p outdir/sample_name --verbose --bamFofn bams.fofn

# Export a bam containing just the reads analyzed for CYP2D6.  
# Reads are labeled by allele and colored by haplotype.
# (use the --grayscale option if you don't like the colors!)
# use logLevel DEBUG for more details on the allele calling
pangu -p outdir/sample_name --verbose -x --logLevel DEBUG  <inBam>  

# Default calling mode is 'wgs'.  
# Use capture/amplicon/consensus modes for other data types
pangu -p outdir/sample_name --verbose -x --mode consensus --logLevel DEBUG  <inBam>
```

## Caller Outputs

Results will be generated with the prefix given, else in the cwd.  Diplotype calls are found in the \<prefix\>[/\_]call.log file.
If you include the `-v,--verbose` option, results will be printed to the screen.
```
pangu -p demo_data/example/NA20129 -x -v demo_data/NA20129_wgs.GRCh38.bam
2022-10-21 18:46:39,921 | INFO | Processing demo_data/NA20129_wgs.GRCh38.bam
2022-10-21 18:46:47,426 | INFO | Identifying SV signatures
2022-10-21 18:46:50,005 | INFO | Phasing and Calling Alleles
2022-10-21 18:47:01,908 | INFO | Calling Diplotype
2022-10-21 18:47:01,908 | INFO | Diplotype called: CYP2D6 *4x2/*5
2022-10-21 18:47:07,317 | WARNING | 4272 bases in variant region with coverage < 3
2022-10-21 18:47:18,207 | INFO | Exporting demo_data/example/NA20129_labeledReads.bam
```
Detailed results can be found in the `<prefix>_report.json` file
```
cat demo_data/example/NA20129_report.json
[
    {
        "input": "demo_data/NA20129_wgs.GRCh38.bam",
        "diplotype": "CYP2D6 *4x2/*5",
        "copynumber": 2,
        "haplotypes": [
            {
                "call": "*4x2",
                "alleles": [
                    {
                        "call": "*4",
                        "num_reads": 24,
                        "meanCover": 16.54,
                        "maxCover": 19,
                        "minCover": 15,
                        "rsIDs": [
                            "rs3892097"
                        ]
                    },
                    {
                        "call": "*4",
                        "num_reads": 13,
                        "meanCover": 11.53,
                        "maxCover": 13,
                        "minCover": 10,
                        "rsIDs": [
                            "rs3892097"
                        ]
                    }
                ]
            },
            {
                "call": "*5",
                "alleles": [
                    {
                        "call": "*5",
                        "num_reads": 10,
                        "meanCover": 1.89,
                        "maxCover": 2,
                        "minCover": 1,
                        "rsIDs": []
                    }
                ]
            }
        ],
        "warnings": [
            "4272 bases in variant region with coverage < 3"
        ]
    }
```
Color:

![demo image](images/NA20129_igv_color.png?raw=true "Labeled Reads")

Grayscale:

![demo image](images/NA20129_igv_grayscale.png?raw=true "Labeled Reads (grayscale)")

## Release History
* 0.1.0 - Initial Release (10/18/2022)
* 0.1.1 - Minor updates to report and amp/capture mode
  * add minFreq for amplicon and capture data variant finding
  * include copy number and warnings in json record
  * add pharmcat-formatted tsv output per sample
  * stricter rules for phasing of duplicates for capture data
* 0.2.0 - Reorganization of src files, minor bug fix
  * Updated install via pip and/or bioconda
  * Fixed copy number counting bug
* 0.2.1 - Include GRCh38 reference in repo for cleaner install


## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

FOR RESEARCH USE ONLY. NOT FOR USE IN DIAGNOSTICS PROCEDURES.
