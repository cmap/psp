dry
===========
Performs filtering and normalization of P100 and GCP data.

Converts level 2 to level 3 data.
Required input is a filepath to a gct file, a filepath to an output file, 
and what to name the output file. Output is writing the gct file.

N.B. This script requires a configuration file. You can specify the location 
of this config file with the optional argument -psp_config_path. Otherwise, 
it will look for psp.cfg in your home directory. See example_psp.cfg for an
example config file.

Example usage:
```
python dry/dry.py input.gct /path/to/output/dir output.gct -sample_nan_thresh 0.9 -subset_normalize_bool -optim -force_assay PR1
```

Maintainer
----------
Lev Litichevskiy	
lev@broadinstitute.org

Additional information/documentation
----------
Specs for dry: https://docs.google.com/a/broadinstitute.com/document/d/1zxKEX_xO_ZrqI4EbogX1by9WRgs0VO24_O15e11Xw58/edit?usp=sharing
