# Script for plotting ts/tv values for QUAL and INFO scores from GATK HaplotypeCaller

This script plots transition to transversion ratios against six INFO stats and QUAL scores emitted by GATK HaplotypeCaller for biallelic SNPs.

## Requirements

The script is written in python and uses [Numpy](https://numpy.org/) and [Matplotlib](https://matplotlib.org/). You can install both of them for example via conda using 

```
conda install conda-forge::matplotlib
conda install conda-forge::numpy
```

## Usage

Basic usage is

```
python3 plot_tstv.py --vcf vcf_input.vcf.gz
```

Where **vcf_input.vcf.gz** is your vcf file, either bgzipped (file ends with **.gz**) or uncompressed (file ends with **.vcf**).

### optional parameters

 + **--n_sites** (deault 1000000) The size of the preallocated arrays, in which all input data will be saved in memory. If this number is smaller than the number of biallelic SNPs in your dataset, the script will have to increase array size by 1000000 once the array is full, which may slow the script down a bit. If it is larger the script will use more memory than necessary.
 + **--plot_n_windows** (default 100) The number of windows into which your sorted data will be split
 + **--output_prefix** String, which will be appended to the beginning of all output files

## Output

The script will produce 8 plots in pdf formats, 6 for the GATK INFO stats (FS.pdf, MQ.pdf, MQRankSum.pdf, QD.pdf, ReadPosRankSum.pdf and SOR.pdf) and 2 for QUAL, one with a regular Y axis (QUAL.pdf) and another one with a logarithmic Y axis scale (QUAL_logY.pdf).

An example output file for MQ could look as follows:

![Example plot for MQ](example/MQ.pdf)

 + The X axis shows the cumulative proportion of the biallelic SNPs after being sorted by the statistic of interest. It is split in 100 windows (by default, can be changed by setting the **--plot_n_windows** parameter)
 + The left Y axis (red lines) shows ts/tv values for each window
 + The right Y axis (blue lines) shows the mean value for each window for the statistic of interest
 + The dotted green line shows the [GATK hard filtering threshold for the INFO score of interest](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants). Note that the threshold may be either above or below that line, depending on the INFO stat


