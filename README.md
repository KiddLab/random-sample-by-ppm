# random-sample-by-ppm
Randomly sample genomic positions based on a Position Probability Matrix.  

By default, randomly choses a genomic strand.  Also includes scripts for calculating a
PPM based on observed motifs.




## Usage ##

First, create a position probability matrix from a set of observed motifs.  By default,
a pseduo count total of one (0.25 per each of 'A','C','G','G') is added to each column

```
python create-ppm.py --in Gilbert-L1-targetsite.PMID16107723.txt \
--out Gilbert-L1.ppm.txt --pseduo 1
```