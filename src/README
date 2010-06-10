Usage:

zcat gt_chrXX.xml.gz | ./gpc <threshold> > output.txt

The <threshold> is optional and specified in percentages. If at least one 
population has a minor allele frequency higher than this threshold, a variant
is reported. Note that this threshold must be met for both nucleotides.

The file gt_chrXX.xml.gz is one of the files in:
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/genotype/
The program is tested with these files, but in principle, any of the files in:
ftp://ftp.ncbi.nih.gov/snp/organisms/*/genotype/
should work.

Make sure the libexpat1-dev library is installed (or any other library that
provides expat.h and expat.lib).

Note that this program uses "37:GRCh37" as a reference sequence, if this is
not correct, change the REF_BUILD constant in the source.

The output is structured as follows:

Name       | type    | Description
---------------------------------------------------------------
position   | integer | Position of the first SNP.
id1        | integer | dbSNP rs number of the first SNP.
id2        | integer | dbSNP rs number of the second SNP.
freq_id2_1 | float   | Frequency of allele 1 of the second SNP.
freq_id2_2 | float   | Frequency of allele 2 of the second SNP.
freq_id1_1 | float   | Frequency of allele 1 of the first SNP.
freq_id2_2 | float   | Frequency of allele 2 of the first SNP.
