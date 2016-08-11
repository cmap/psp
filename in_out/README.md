GCToo
-----
GCToo is an implementation in Python of storing gct files as objects in memory.
GCToo objects are built on [pandas](http://pandas.pydata.org/pandas-docs/version/0.17.0/index.html#) dataframes.  

This directory also contains concat_gctoo.py, which can be used for concatenating gct files together.

The files in this directory are copies of files in the cmap/io/GCToo directory in the 
[pestle repository](https://github.com/cmap/pestle) maintained by the Connectivity Map.

Example of using concat_gctoo
-----------------------------

1) You have a bunch of files that start with 'LINCS_GCP' in /some/dir that
you want to concatenate. From the top directory of this repo, run the
following command in your Terminal:

```
python in_out/concat_gctoo.py /some/dir/LINCS_GCP* concated.gct
```

2) You have 2 GCToo objects in memory that you want to concatenate.

```
import in_out.concat_gctoo as cg
concated = cg.hstack([gct1, gct2])
```


Maintainer
-------
Lev Litichevskiy  
lev@broadinstitute.org  
August 2016
