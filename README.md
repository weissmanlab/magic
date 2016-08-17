MAGIC
=======

This is the implementation for Minimal Assumption Genomic Inference of Coalescence (MAGIC).
From genome sequences, MAGIC can be used to infer the distribution of pairwise coalescence times
(i.e., "_N_<sub>e</sub>(_t_)"), as well as the distribution of lengths of any part of the coalescent trees
(e.g., total branch length, tips, branches above cherries, etc).
MAGIC works with the [msmc-tools](http://github.com/stschiff/msmc-tools) set of programs to process the initial data.
Everything, including the parts of `msmc-tools` that MAGIC interacts with, is written
in Python, and is designed to be modular, so you are encouraged to edit and add to it.
(Like `msmc-tools`, MAGIC requires Python 3; it also requires [SciPy](https://www.scipy.org).)
Here's the basic workflow:

## combo\_prep.py

This is an enhanced version of `msmc-tools`' [generate\_multihetsep.py](https://github.com/stschiff/msmc-tools/blob/master/generate_multihetsep.py).
It takes the same input formats and returns the same output, ready for MSMC,
and is used with *almost* the same syntax (see the [msmc-tools user guide](https://github.com/stschiff/msmc-tools/blob/master/README.md)).
The only substantive difference is that if the optional keyword argument `--coverfile <filename>` is used,
it additionally returns a file giving the sequence coverage of the whole genome in short windows,
which MAGIC needs.
The other slight difference is that the `--mask` and `--negative_mask` arguments have been replaced
by `--masks` and `--negative_masks`, respectively, which work slightly differently.
Whereas with `generate_multihetsep.py` you would enter multiple masks as `--mask mask1 --mask mask2 ...`,
with `combo_prep.py` this would be `--masks mask1 mask2 ...`, which reduces redundancy and works better with globbing, etc.
Note that this means that these arguments must come *after* the vcf files in the command line -- see the example in the [Summary](#summary) below.

## windower.py

This performs the first real step of MAGIC's analysis: it converts the SNP data into a list of
diversity histograms at a range of genomic length scales. 
Specifically, for each length scale (eg, 1kb), it splits the genome into windows of that length
and returns how many of those windows contain 0 SNPs, how many contain 1, etc.
This is the stage where you specify which coalescence time distribution you're trying to find;
`windower.py` will filter all the SNPs and only look at the ones that are relevant for that time.
So, for instance, if you want to find the distribution of the total length of all branches above cherries,
it would count all doubleton SNPs, while if you wanted the distribution of the sum of the tips above 
individual 1's two haplotypes, it would count all of individual 1's singleton SNPs.
There are a few default options built in, but the idea is that you can easily write your own.
The defaults all work with unphased, unpolarized data, so if you want to take advantage of the 
extra information provided by phasing or polarization, you'll have to do it!

The basic syntax is:

	./windower.py <SNP file> --stat <statistic to calculate>

## magic.py

This performs the remaining steps of the analysis. The basic syntax is:

	./magic.py <histogram files> --out <output_prefix>
	
where the histogram files are the outputs of `windower.py`. 
This will produce a file with the suffix `_final.txt` with the best fit
mixture of [gamma distributions](https://en.wikipedia.org/wiki/Gamma_distribution)
for the desired coalescence time distribution. The first column of the file gives
the weights of each component gamma distribution, the second column gives the shape
parameters, and the third gives the scale parameters. 
If the last row has only a single column, it is the estimated probability that the 
desired time is exactly 0. 
(This is expected for many tree features, e.g.,
branches that lie over exactly 1/3 of the haplotypes in the sample; 
for these features you can force `magic.py` to include this line using the optional argument `--zero`.)
Additionally, `magic.py` produces a file with the suffix `_full.txt` that gives more information
about how good the fit is, and a file with the suffix `_LT.txt`
giving the Laplace transform of the distribution at a discrete set of points.
These points are the underlying data used to fit the mixture of gamma distributions, 
and are useful in their own right for model checking.

## magicplots.ipynb

If you have [Jupyter](http://jupyter.org/) installed, this is a template for making a notebook
with plots of the distributions produced by MAGIC. 
It's also useful for doing further quantitative analysis of the results.

# Summary

Here's a basic example of how the remaining steps might look. 
Say that you have used `msmc-tools` to prepare zipped vcf files of chromosome 1 from a sample of many individuals,
with each file having a name like `<i>_chr1.vcf.gz`, 
along with masks for covered sites with names like `covered_<i>_chr1.bed.txt.gz`,
where `<i>` is the label for each individual.
Now you want to analyze this sample.
You would run:

	./combo_prep.py *_chr1.vcf.gz --coverfile chr1_cover.txt --masks covered_*_chr1.bed.txt.gz > chr1.txt
	
Now you can calculate all the distributions you want without re-running this step.
To calculate the total branch length distribution (the default), you would run:

	./windower.py chr1
	./magic.py chr1_counts.txt --out chr1
	
If you've processed multiple chromosomes with `combo_prep.py` and `windower.py` and want to 
combine them all to find the genomic distribution, you would run:

	./magic.py chr*_counts.txt --out genome
	
If you then want to find the distribution of, say, total tip lengths, you would run:

	./windower.py chr1 --stat 1
	./magic.py chr1_1_counts.txt --out chr1_ttip
	
See the code for `windower.py` for the names of other built-in statistics.


	