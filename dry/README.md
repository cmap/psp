dry
===========
Performs filtering and normalization of P100 and GCP data.

Converts level 2 to level 3 data. Required input is a filepath to a gct file
and the path to where output should be saved. The two outputs are a gct file
and a pw file ("plate-well" file, contains QC results).

N.B. This script requires a configuration file. You can specify the location 
of this config file with the optional argument -psp_config_path. Otherwise, 
it will look for psp_production.cfg in your current directory.

Example usage:
```
python dry/dry.py input.gct /path/to/output/dir -out_name output.gct -out_pw_name output.pw -sample_nan_thresh 0.5
```

Maintainer
----------
Lev Litichevskiy	
lev@broadinstitute.org

Additional information/documentation
----------
Specs for dry: https://docs.google.com/a/broadinstitute.com/document/d/1zxKEX_xO_ZrqI4EbogX1by9WRgs0VO24_O15e11Xw58/edit?usp=sharing
