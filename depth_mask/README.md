# Create depth mask for VCF files
This logic behing this script is similar to this approach from [Marek Slenker's GitHub](https://github.com/MarekSlenker/vcf-mask), but it is faster and more straightforward to use.
Basically, it takes a VCF file as input, determines the mean and standard deviation of sequencing depth (DP in the genotype field of the VCF file). It will then create a mask for sites where a pre-defined (via -n or --nsamples) number of samples has depth greater than mean+2*Std.
