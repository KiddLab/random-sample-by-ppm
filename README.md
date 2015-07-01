# random-sample-by-ppm
Randomly sample genomic positions based on a Position Probability Matrix.  

By default, randomly choses a genomic strand.  Also includes scripts for calculating a
PPM based on observed motifs.

Assumes that motifs and PPM correspond to endonuclease cleavage sites reported as in
Gilbert et al. 2005  (http://www.ncbi.nlm.nih.gov/pubmed/16107723).

## Usage ##

First, create a position probability matrix from a set of observed motifs.  By default,
a pseduo count total of one (0.25 per each of 'A','C','G','G') is added to each column

```
python create-ppm.py --in Gilbert-L1-targetsite.PMID16107723.txt \
--out Gilbert-L1.ppm.txt --pseduo 1
```


Next, create the desired number of random sample sets (--n_sets) each consisting of the
desired number of random samples (--n_per_set).  Note, that the sequence of the genome
is read into memory, so be sure that your computer has enough RAM available.  First three
columns of the output will be as a bed file.

```
python sample-based-on-ppm.py \
--ppm Gilbert-L1.ppm.txt \
--genome genomes/hg19.fa \
--n_per_set 100 \
--n_sets 5 \
--outpre sample100.test
```

Optionally, include a bed file of regions to exclude placement (such as near other repeat elements)
using --exclusionbed



