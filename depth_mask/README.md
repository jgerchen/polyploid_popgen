# Create depth mask for VCF files

## Principle

The logic behing this script is similar to this approach from [Marek Slenker's GitHub](https://github.com/MarekSlenker/vcf-mask), but it is faster and more straightforward to use.
Basically, it takes a VCF file as input, determines the mean and standard deviation of sequencing depth (DP in the genotype field of the VCF file). It will then create a mask for sites where a pre-defined (via -n or --nsamples) number of samples has depth greater than mean+2*Std.

## Requirements

### Software
Python3 with NumPy

### Input files
A vcf file, which can be gzipped (determined if the ending of the file name is .gz), with snps and invariants both after GATK filtering (filtered, but not removed), but without insertions, deletions and spanning deletion sites (i.e. the sites with *)
