# Filtering subsets of VCF files

## The issue

Since we usually want to work with biallelic SNPs, we typically filter out multi-allelic sites at some point of our variant calling pipelines. However, the number of multi-allelic sites increases with larger datasets, because there will be additional singleton overlapping biallelic sites, which will then be lost for further analyses. The severity of this issue depends on the size of datasets as well as several demographic parameters of the included populations.

## Potential solutions

For many analyses we are fine with looking at just a subset of our complete dataset at a time. So for example if we want to look at population-level stats like pi or Tajima's D we may look at a single population at a time. If we want to do estimate pairwise stats like FST we could look at subsets of a specific population pair. Generating subsets of our datasets so that we include only the individuals we are interested in could also help recover biallelic sites that would otherwise be lost, because our subsets will not include all alleles for many multi-allelic sites and we could recover multi-allelic sites for which we have only two alleles in our subset.

However, while this seems straightforward there are some issues with using and filtering these subsets. We usually filter sites based on several properties on both the genotype and the site level. One important set of information for filtering are INFO stats emitted by our variant callers (e.g. GATK HaplotypeCaller or freebayes). Importantly, INFO stats are based on the whole dataset that went into the variant calling pipeline and give us a set of statistics that are supposed to help us determine if called variation at a specific site is likely to be caused by a true variant anywhere in the dataset. Generally, the more data goes into the variant caller, the more reliable INFO stats will be. Also all but the most trivial INFO stats cannot be recalculated without rerunning the variant caller.

Another important issue with INFO stats is that stats for SNPs are not directly comparable with INFO stats for indels or complex variants. And many of the multi-allelic sites that we want to recover are going to be some form of complex variant, were some rare indel will overlap with a SNP. If we just ignore the indel and look at the SNP our INFO stat will not reflect information on the SNP, but on the complex variant as a whole.

### Rerunning the variant calling based on gvcf files

If you are interested in multiple subsets of your data you could just do the variant calling separately for each of those subsets. For GATK HaplotypeCaller you could run the first steps (aligning and filtering reads, running the first steps on GATK HaplotypeCaller) on your complete dataset and then do the actual variant calling (GenomicsDBImport and GenotypeGVCFs) for your subset. Then you would run filtering on ech of those subsets.  

**Pros**:

- Our dataset will only include variants found in our subset and may automatically recover many bi-allelic SNPs, that would otherwise be marked as multi-allelic.

- INFO stats are properly recalculated for the subset in question

**Cons**:

- You have to redo the variant calling for every subset! This can be very computationally intensive.

- INFO stats will only reflect the variation in the subset and information on the whole dataset will be lost. As a result you may end up filtering some low-frequency variants (in your subset) that you may not filter otherwise if you had information on the full dataset.

### Filtering the raw VCF file based on all your data

Another option when working with subsets of data is to use the raw unfiltered VCF file based on your complete dataset and then subset it and filter data for your subsets. Generally all of this can be done using bcftools, but there are some fixes we have to apply when working with mixed-ploidy datasets. I describe the basic steps and bcftools commands below:

We start by subsetting our raw VCF file using bcftools


```
bcftools view -S subset_samples.list -o raw_vcf_subset.vcf.gz -O z raw_vcf.vcf.gz
```

Here **subset_samples.list** is a simple text file containing the names of the samples (must be the same as in the VCF header) that we want to keep for our subset. An example can be found [here](alnus_samples.list).

If we work with polyploids, we have to fix the header of our VCF file or bcftools will crash at the following steps.

First we extract the VCF header

```
bcftools view -h raw_vcf_subset.vcf.gz > raw_vcf_subset.header
```

Then we use sed to replace some text in the header

```
sed -i s/ID=PL,Number=G/ID=PL,Number=./g raw_vcf_subset.header
```

Then we use bcftools to make a new VCF file with the updated header (called raw_vcf_subset.**rh**.vcf.gz).

```
bcftools reheader -h raw_vcf_subset.header -o raw_vcf_subset.rh.vcf.gz raw_vcf_subset.vcf.gz 
```

Now we fixed our VCF file for the next steps. If we want to filter individual genotypes by depth now would be a good point to do that.

```
bcftools filter -O z -o raw_vcf_subset.rh.gtdepth.vcf.gz  -i 'FMT/DP>4' --set-GTs . raw_vcf_subset.rh.vcf.gz 
```

This will replace all genotypes that are covered by four reads or less to no-call and create a new VCF file called raw_vcf_subset.rh.**gtdepth**.vcf.gz .

Now we can run the magic bcftools command that will look at the REF and ALT columns in our VCF file and only retain alleles that are found in our dataset. So it will remove all alleles that are not found in our subset, including those that got removed by depth filtering of genotypes in the previous step.

```
bcftools view -a -o raw_vcf_subset.rh.gtdepth.sub.vcf.gz raw_vcf_subset.rh.gtdepth.vcf.gz
```

And now we have our properly subset vcf file on which we can run further filtering steps. Importantly, i**the INFO stats in your VCF file will still be the same as in the raw VCF file for the whole dataset**. This has some consequences for filtering your data:

- For SNPs that were biallelic in the complete dataset this may be actually desirable, because your INFO stats will be based on a greater set of information than your subset alone.

- However, for recovered SNPs that are now biallelic the INFO stats will not be useful, since they are based on the multi-allelic information.

**Pros**:

- You don't have to redo the variant calling!

- Your INFO stats will reflect information from your complete dataset for bi-allelic sites.

**Cons**:

- Your INFO stats will not accurately reflect information for your recovered sites (which were the reason why you did the subsetting in the first place)

## How to filter

As you can see we are in a bit of a dilemma here: if we redo the variant calling for subsets of data our INFO stats may be less informative in general, while if we subset our raw VCF our INFO stats will be accurate for biallelic SNPs, but maybe not useful for the multi-allelic sites we wanted to recover in the first place. So what can we do?

- **Ignore the issue and just do the filtering for the subset as we would do for the complete dataset**

This would be the obvious way to go if you recalled your VCF for subsets of data. If you subset your raw VCF the question would be to what degree using INFO stats from the multiallelic sites for diploids would bias your filtering. After all we don't gain anything if we recover variants by subsetting and then discard the recovered variants again in the follow up filtering, because the INFO stats based on multi-allelic sites don't make the filtering thresholds. However, I don't know how big this issue would be so this **would have to be tested**.

- **Not use INFO stats for filtering**

While using the [GATK best practices filtering thresholds](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) are a kind of standard commonly applied to non-model species, this is far from ideal. These threshold were developed for human data, and besides a number of other shortcomings humans are boring diploids. I think that for most species these thresholds may be overly conservative and remove quite a lot of sites, which may be true positives. This may be alleviated to some degree by [looking at ts/tv plots and adjusting thresholds manually](https://github.com/jgerchen/polyploid_popgen/tree/main/tstv_plot), which may improve recovery of sites substantially, or if training data are available machine learning approaches like [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR) could be used.

However, all of this is not useful for subsets of data and the question is if you would introduce too strong biases if you did not use INFO scores for filtering at all. I think for many analyses this may actually be the better choice. After all most false positives are going to be due to sequencing errors, which will be at low frequencies, which in turn are expected to have little effect on statistics like FST. Furthermore, many of these errors will be at low coverage sites which will already be removed by depth-filtering genotypes. In the end the question would be how the higher number of recovered sites would balance against the introduced number of false positives, which is very difficult to answer empirically.

- **Filter only sites that were biallelic in the complete dataset based on INFO scores and do something else for the recovered sites**

This may sound like the goldilocks option, but it would also imply that the two classes of sites are treated differently, which may not be ideal. Also it would require some additional steps for subsetting datasets, but this should be doable with some fiddling with bcftools.
