# Find variable CpG loci

## Installation
nstall the [expat](http://expat.sourceforge.net/) development files:

    apt-get install libexpat1-dev

Compile the program:

    make

## Usage
To run the program:

- Get one of the files from
  [this](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/genotype/)
  location.
- If we assume that the downloaded file is named `gt_chrXX.xml.gz`, use
  the following command to find the CpG candidates:

    zcat gt_chrXX.xml.gz | ./gpc <threshold> > output.txt

The `threshold` is optional and specified in percentages. If at least one 
population has a minor allele frequency higher than this threshold, a variant
is reported. Note that this threshold must be met for both nucleotides.

Although we have not tested this program on files related to other species, we
see no reason why it should not work. To find input files for other species, go
[here](ftp://ftp.ncbi.nih.gov/snp/organisms/), select one of the directories
and choose the genotype subdirectory (if present).

Note that this program uses `37:GRCh37` as a reference sequence, if this is
not correct, change the `REF_BUILD` constant in the source.

The output is structured as follows:

|Name       |type    |Description
|---        |---      |---
|position   |integer |Position of the first SNP.
|id1        |integer |dbSNP rs number of the first SNP.
|id2        |integer |dbSNP rs number of the second SNP.
|freq_id2_1 |float   |Frequency of allele 1 of the second SNP.
|freq_id2_2 |float   |Frequency of allele 2 of the second SNP.
|freq_id1_1 |float   |Frequency of allele 1 of the first SNP.
|freq_id2_2 |float   |Frequency of allele 2 of the first SNP.
