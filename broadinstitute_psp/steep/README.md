steep
===========
Computes similarity (e.g. Spearman correlation).

If given 1 gct, steep will compute all pairwise similarities between its columns.
If given 2 gcts, steep will compute pairwise similarities between the columns
of gct1 and the columns of gct2, but not the similarities within a single gct.
Required arguments are 1 gct, the output directory, and the name for the
output similarity gct.

Example usage:
```
python steep/steep.py -i /path/to/first.gct -i2 /path/to/second.gct
```

Maintainer
----------
Lev Litichevskiy	
lev@broadinstitute.org  
January 2017

