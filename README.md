# simple-fst
Perl scripts to calculate per snp FST

## what?

These are some simple scripts written in perl, from early in my time at MPI/EVA leipzig.

pairwise_FST_WC.pl takes a vcf, also in gz format, and calculates per snp fst,
as defined by Weir + Cockerham 1984 (Evolution 38(6) pp1358-1370, 1984.). 
Two txt files, containing the sample IDs, one each for each population.

The main reason I wrote these instead of using vcftools, was so that I could also return the estimated a, b, and c parameters
for more flexible down stream analyses e.g. windowed average fst.

pairwise_FST_patterson.pl does the same thing, more or less, but instead uses the fst derivation in Bhatia et al Genome Res. 2013. 23: 1514-1521.
It returns the numerator (nHat) and denominator (dHat), as well as the fst.

Note: I tested pairwise_FST_WC.pl against output from vcftools, and they were identical, but of course please confirm for yourself, using the small test data set.

Further note: These are fairly basic, and if i had more time I would rewrite them but they do the job as is.

## how do?

perl scripts/pairwise_FST_WC.pl data/kGenomes_chr14_small.vcf.gz data/MSL_61.txt data/MXL_61.txt 14 results

```there might be some perl module dependnecy issues depending on your local perl installation!```
