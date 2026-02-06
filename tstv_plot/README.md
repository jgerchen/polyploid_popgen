# Script for plotting ts/tv values for QUAL and INFO scores from GATK HaplotypeCaller

This script plots transition to transversion ratios against six INFO stats and QUAL scores emitted by GATK HaplotypeCaller for biallelic SNPs.

## Requirements

The script is written in python and uses [Numpy](https://numpy.org/) and [Matplotlib](https://matplotlib.org/). You can install it for example via conda using 

´´´
conda install conda-forge::matplotlib
conda install conda-forge::numpy
´´´

## Usage

Basic usage is

´´´
python3 plot_tstv.py --vcf vcf_input.vcf.gz
´´´

Where **vcf_input.vcf.gz** is your vcf file, either bgzipped (file ends with **.gz**) or uncompressed (file ends with **.vcf**).

### optional parameters

 + **--n_sites** (deault 1000000) The size of the preallocated arrays, in which all input data will be saved in memory. If this number is smaller than the number of biallelic SNPs in your dataset, the script will have to increase array size by 1000000 once the array is full, which may slow the script down a bit. If it is larger the script will use more memory than necessary.
 + **--plot_n_windows** (default 100) The number of windows into which your sorted data will be split
 + **--output_prefix** String, which will be appended to the beginning of all output files

## Output
