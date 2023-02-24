# Create depth mask for VCF files

## Principle

The logic behing this script is similar to this approach from [Marek Slenker's GitHub](https://github.com/MarekSlenker/vcf-mask), but it is faster and more straightforward to use.
Basically, it takes a VCF file as input, determines the mean and standard deviation of sequencing depth (DP in the genotype field of the VCF file). It will then create a mask for sites where a pre-defined (via -n or --nsamples) number of samples has depth greater than mean+2*Std (calculated for each sample).

## Requirements

### Software
Python3 with NumPy

### Input files
A vcf file, which can be gzipped (determined if the ending of the file name is .gz), with snps and invariants both after GATK filtering (filtered, but not removed), but without insertions, deletions and spanning deletion sites (i.e. the sites with *)

## Usage

### Parameters

**-v, --vcf**

Input VCF file (gzipped or not)

**-o, --output**

Output file containing sites to be masked (contig and position in columns, tab separated)

**-c, --countout**

Output file containing all sites and the number of samples that are above cutoff for filtering (contig, position and number of samples in columns, tab separated)

**-p, --plothist**

Output file containing histogram with summed-up numbers of counts of samples that are above cutoff for all sites (number of samples and counts in columns, tab separated)

**-n, --nsamples**

Minimum number of samples above cutoff for a site to be masked

**-l, --locinumber**

Number of loci in the VCF file. Ideally this should be the exact number from the VCF file. If the number is higher, the program will need more memory, if it is lower it will have to reallocate memory to the internal genotype array, which is computationally inefficient.
For a gzipped VCF file, this can be determined using 
```
zcat vcf_file.vcf.gz | grep -c "^[^#]"
```
For an unzipped VCF file it would be 

```
grep -c "^[^#]" vcf_file.vcf
```
### Example usage

The program could be run like
```
python3 make_depth_mask.py -v vcf_file.vcf.gz -o out_file.tsv -c counts_out.tsv -p hist_out.tsv -n 18 -l 69245
```

Based on the histogram in hist_out.tsv (determined by -p), you may decide to choose a different value for -n. For this it is not necessary to rerun the script, but you can just filter the counts_out.tsv file using the following code (here for n>=15)
```
awk '{if ($3>=15) print $1"\t"$2}'
```
### Transforming output files

The output files are a simple list of single genomic positions. You can transfer this into bed format using
```
awk '{printf "%s\t%d\t%d\n" ,$1,($2 - 1),$2}' out_file.tsv > output.bed
```
Then you can use bedops to merge adjacent variants into continuous stetches (which will be preferable for use with e.g. GATK)
```
bedops -m output.bed > output_merged.bed
```


