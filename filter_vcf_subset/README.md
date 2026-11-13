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

















Pros:

You don't have to redo the variant calling!

Cons:
