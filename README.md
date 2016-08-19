MAGIC
=======

This is the implementation for Minimal Assumption Genomic Inference of Coalescence (MAGIC).
From genome sequences, MAGIC can be used to infer the distribution of pairwise coalescence times
(i.e., "_N_<sub>e</sub>(_t_)"), as well as the distribution of lengths of any part of the coalescent trees
(e.g., total branch length, tips, branches above cherries, etc).
MAGIC works with the [msmc-tools](http://github.com/stschiff/msmc-tools) set of programs to process the initial data.
Everything, including the parts of msmc-tools that MAGIC interacts with, is written
in Python, and is designed to be modular, so you are encouraged to edit and add to it.

## Requirements
- [Python 3.4+](http://www.python.org)
- [Scipy](http://www.scipy.org/) (including Numpy)
- [msmc-tools](http://github.com/stschiff/msmc-tools)
- Optional, just for plotting and interactive analysis:
	- [Jupyter](http://www.jupyter.org)
	- [matplotlib](http://matplotlib.org/)

[Anaconda](http://www.continuum.io) is an easy way to install Python, Scipy, Jupyter, and matplotlib together, along with a lot of other stuff.


## Basic workflow

1. Pre-process the data using [msmc-tools](http://github.com/stschiff/msmc-tools) to produce vcf files and masks.
2. Use combo_prep.py to convert these into a list of SNPs and a sequencing coverage file.
3. Use [windower.py](#windower.py) to filter the SNPs for those that are informative about the desired coalescence times, and find their distribution across genomic windows.
4. Use [magic.py](#magic.py) to infer the distribution of coalescence times from the window diversity distributions.
5. Plot the results with [magicplots.ipynb](#magicplots.ipynb).

See the [example](#example) for specifics.

## Example

Here's a basic example of how you might use MAGIC. 
Say that you have used `msmc-tools` to prepare zipped vcf files of chromosome 1 from a sample of many individuals,
with each file having a name like `<i>_chr1.vcf.gz`, 
along with masks for covered sites with names like `covered_<i>_chr1.bed.gz`,
where `<i>` is the label for each individual.
(Not sure how to do this? [Instructions.](https://github.com/stschiff/msmc-tools/blob/master/README.md) Still not sure? [Ask!](https://groups.google.com/forum/#!forum/msmc-popgen))
Now you want to analyze this sample.
You would run:

	./combo_prep.py *_chr1.vcf.gz --coverfile chr1_cover.txt --masks covered_*_chr1.bed.gz > chr1.txt
	
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
	
[See the other built-in statistics that you can estimate.](#specifying-features-to-estimate)

To look at all the distributions, run:

	jupyter notebook
	
and open `magicplots.ipynb`.
Change the variable names to `stats = ["chr1", "genome", "chr1_ttip"]` and `pairstats = ["chr1", "genome"]` and run all the cells.





# combo\_prep.py

This is an enhanced version of msmc-tools' [generate\_multihetsep.py](https://github.com/stschiff/msmc-tools/blob/master/generate_multihetsep.py).
It takes the same input formats and returns the [same list of SNPs](https://github.com/stschiff/msmc/blob/master/guide.md#input-file-format), ready for MSMC.
See the [msmc-tools user guide](https://github.com/stschiff/msmc-tools/blob/master/README.md) for instructions
-- everything is exactly the same, with two exceptions.
The substantive one is that if the optional keyword argument `--coverfile <filename>` is used,
`combo_prep.py` returns an additional file giving the sequence coverage of the whole genome in short windows,
which MAGIC needs.
The other difference is superficial: the `--mask` and `--negative_mask` arguments have been replaced
by `--masks` and `--negative_masks`, respectively, which work slightly differently.
Whereas with `generate_multihetsep.py` you would enter multiple masks as `--mask mask1 --mask mask2 ...`,
with `combo_prep.py` this would be `--masks mask1 mask2 ...`, which reduces redundancy and works better with globbing, etc.
Note that this means that these arguments must come *after* the vcf files in the command line -- see the [example](#example).

# windower.py

This performs the feature-specific part of the data preparation. 
Given a desired coalescence time feature (e.g., pairwise time or total branch length), it first filters all the SNPs for the ones that are relevant for that time.
Then, for each of a range of length scales, 
it splits the genome into windows of that length
and returns how many of those windows contain 0 SNPs, how many contain 1, etc.

The basic syntax is:

	./windower.py <SNP file> --stat <statistic to calculate>
	
This produces a file (by default with suffix `_counts.txt`) with one line for each length scale, starting from the shortest.
Each line consists of a number of comma-separated entries $$$n\ k_n$$$, where $$$n$$$ is a number of SNPs and $$$k_n$$$ is the number of windows of that length that contain $$$n$$$ SNPs.

### Specifying features to estimate

For each piece of a potential coalescent tree, mutations that occur on that piece
will show up as SNPs with a characteristic distribution across individuals.
For example, the set of all branches above cherries corresponds to doubleton SNPs, while the tips above individual 1's two haplotypes correspond to individual 1's singleton SNPs.

There are a few default options for tree features built into `windower.py`:

- `tbl` (default): Total branch length (all polymorphisms). To find the pairwise coalescence time distribution, run `windower.py` on each (diploid) individual with this option, and then run `magic.py` on all the outputs together. (Note the factor of 2 difference between the pairwise T<sub>MRCA</sub> and the pairwise TBL though.)
- `indiv_tips`: Length of tips above each individual (singletons, separated by individual). For a sample of $$$n$$$ individuals, this will return $$$n$$$ output files, one for each individual. Run `magic.py` separately on each of these outputs to find the distribution of tip lengths above each individual, or on all of them together to find the combined distribution across all individuals.
- Integer $$$j$$$: Length of branches lying above exactly $$$j$$$ haplotypes in the unrooted tree (SNPs where minor allele has $$$j$$$ copies). Running `magic.py` on the output will give a more powerful version of the (folded) site-frequency spectrum. While the site-frequency spectrum just gives the mean length of the corresponding set of tree branches, this will give the full distribution of lengths.

The idea, though, is that you can easily write your own filter for the features you want to know about.
The defaults all work with unphased, unpolarized data, so if you want to take advantage of the 
extra information provided by phasing or polarization, you'll *have* to do write your own!



# magic.py

This is where MAGIC infers the coalescence time distributions from the genomic diversity.
The basic syntax is:

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
All times are normalized with the per-base mutation rate, so, for example, $$$0.001$$$ means $$$0.001/\mu$$$.

Additionally, `magic.py` produces a file with the suffix `_LT.txt`
giving the Laplace transform of the distribution at a discrete set of points.
These points are the underlying data used to fit the mixture of gamma distributions, 
and are useful in their own right for model checking.
Each row of the file is of the form 
$$$s \ \ \tilde p(s) \ \ \sigma(\tilde p(s))$$$: a value of the Laplace transform variable, the estimated value of the Laplace transform at that point, and the error in the estimation.

# magicplots.ipynb

If you have [Jupyter](http://jupyter.org/) installed, this is a template for making a notebook
with plots of the distributions produced by MAGIC. 
It's also useful for doing further quantitative analysis of the results.



	