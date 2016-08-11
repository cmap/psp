steep
===========
Computes similarity (e.g. Spearman correlation).

If given 1 gct, steep will compute all pairwise similarities between its columns.
If given 2 gcts, steep will compute pairwise similarities between the columns
of gct1 and the columns of gct2, but not the similarities within a single gct.
Required arguments are 1 gct, the output directory, and the name for the
output similarity gct.

N.B. This script requires a configuration file. You can specify the location
of this config file with the optional argument -psp_config_path.
Otherwise, it will look for psp_production.cfg in your current directory.

N.B. This directory is being refactored, so please ignore old_steep.py for now.

Example usage:
```
python steep/steep.py new.gct /path/to/output/dir/ new_v_old.gct -in_gct2 old.gct
```

Maintainer
----------
Lev Litichevskiy	
lev@broadinstitute.org  
August 2016

More info
---------
Please go to the [specs for steep](https://docs.google.com/a/broadinstitute.com/document/d/10--axSFOdRAKgjSl0zIPhMcQs9bSVyEbHVPgcqCu240/edit?usp=sharing).
TODO(lev): update this spec
