
# check_mendel.pl
Perl script to check Mendelian errors in VCF file using pedigree information

## Usage
`$ perl check_mendel.pl --vcf [VCF] --ped [PED] --out [output_prefix] `

## Input
* VCF file with GT field
* PED file with trio or duo information

## Output
* [output_prefix].tbl : A table with Mendelian error counts for all possible parent/child genotype combinations
* [output_prefix].AC : Number of errors and variants for each allele count bin
* [output_prefix].var : List of all variants with number of Mendelian errors annotated
