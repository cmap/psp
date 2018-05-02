# Proteomics Signature Pipeline (PSP)

PSP is a collection of Python scripts that enable data processing and analysis of GCP and P100 proteomic data produced by the LINCS Proteomic Characterization Center for Signaling and Epigenetics (PCCSE) at the Broad Institute. You can download the raw data from the [Panorama Repository](https://panoramaweb.org/labkey/project/LINCS/begin.view? "Panorama Repository"). (You will want the unprocessed gcts.)

In addition to the scripts available here, you can get quick access to these data with a suite of web apps. Check them out at [clue.io/proteomics](https://clue.io/proteomics "Proteomics on Clue").

|Production | [![Build Status](https://jenkins.clue.io/buildStatus/icon?job=TEST-PSP-MASTER)](https://jenkins.clue.io/job/TEST-PSP-MASTER) |
--- | --- |
|Dev | [![Build Status](https://jenkins.clue.io/buildStatus/icon?job=TEST-PSP)](https://jenkins.clue.io/job/TEST-PSP) |


## Maintainer

Lev Litichevskiy  
lev@broadinstitute.org  
May 2018

## Setting up your environment

### To set up your environment the first time:

  1. First, you must clone this repository into a local directory on your computer. For example, I cloned the repository in the directory `/Users/lev/code/psp`. If you need more information about cloning, go to this page provided by Github: https://help.github.com/articles/cloning-a-repository/.

  2. To manage our Python environment, we'll use a program called conda. Download conda from the following website: http://conda.pydata.org/miniconda.html. Miniconda or Anaconda will do, but we recommend Miniconda because it's more lightweight.

  3. Now, we will continue with the setup using the command line, so open the terminal or command prompt.

  4. Type `conda info` to verify that conda has been installed on your system. You should see some information about the "Current conda install." If you don't, then conda has not been installed properly.

  5. We will now create an environment with conda that will allow us to use PSP. If on OSX or Linux, type the following in your Terminal:

      ```
      conda create --name psp_env --channel bioconda python=2 pandas scipy h5py cmappy=3.2.0
      ```
      
      'psp_env' will be the name of your conda environment, and the things after it are the packages that our environment will contain. Note that we are using python2, not python3. We also have to specify that we should also look in the `bioconda` channel in order to find the [cmapPy](https://github.com/cmap/cmappy "cmapPy Github") package (tools for interacting with .gct and .gctx files), namely the 3.2.0 version. You'll have to type 'yes' to proceed through the installation.
      
      If on PC, type the following in your terminal:
      
      ```
      conda create --name psp_env --channel bioconda python=2 pandas scipy h5py
      activate psp_env
      pip install cmappy==3.2.0
      ```
      
      Unfortunately, bioconda (which is where [cmapPy](https://github.com/cmap/cmappy "cmapPy Github") is hosted) does not support Windows, so we have to use pip to install it.
      
  6. OPTIONAL: There are additional dependencies that you will need in order to use the tasseography scripts (i.e. showing connections as graphs):
    
      ```
      conda install --channel conda-forge python-igraph matplotlib
      ```
      
      igraph is a package for manipulating graphs, and you'll need matplotlib to produce output figures.

  7. To activate your environment, type `source activate psp_env`, or if you are on a PC, `activate psp_env`. You should now see `[psp_env]` or `(psp_env)` prepended to the start of your command prompt. For example:

      ```
      (psp_env) /Users/lev/code/psp $
      ```

  8. We will now run a script to make our environment aware of the contents of the PSP repository that we cloned. Move to the top directory of your cloned repository (e.g. `/Users/lev/code/psp`), and then type:

      ```
      python setup_psp.py develop
      ```   
  
  9. That's it! To make sure that everything has been set up correctly, navigate to the `broadinstitute_psp` directory of PSP and try executing one of the Python test scripts:

      ```
      cd /Users/lev/code/psp/broadinstitute_psp/
      python dry/test_dry.py
      ```
      
      To be thorough, you can run all tests with the following command:
      
      ```
      python -m unittest discover ./ -p "test_*.py"
      ```
      
      Test files should be executed from the `broadinstitute_psp` directory, but other scripts can be run from anywhere. For example, you could run dry from anywhere by supplying the full path to the script. For example:
      
      ```
      python /Users/lev/code/psp/broadinstitute_psp/dry/dry.py --help
      ```
      
      However, before doing this, make sure to check out the section below about the configuration file.
  
### To set up your environment after the first time:

  1. Activate your conda environment by typing `source activate psp_env`, or if you are on a Windows computer, `activate psp_env`.
  2. It's easiest to run your code from the `broadinstitute_psp` directory (the one that contains the dry, steep, etc. directories), so move there and try executing one of the test scripts again:
    
      ```
      cd /Users/lev/code/psp/broadinstitute_psp
      python dry/test_dry.py
      ```

## cmapPy

cmapPy is the repository of Python tools for interacting with .gct and .gctx files. This repo relies heavily on it. We install it using conda or pip above. Note that you also get several command line tools for free, such as gct2gctx, gctx2gct, subset, and concat. For example, type the following in your terminal:

`concat -h` 

If you see the help page for concat, you can use this tool directly from the command line. See the [cmapPy repo](https://github.com/cmap/cmappy "cmapPy Github")  for more information.


## Configuration file

Several of the files in PSP -- notably dry.py -- require a configuration file. This file contains various parameters used by the production code (e.g. what values to consider NaN, default thresholds, provenance code abbreviations, where to find certain files, etc.). By default, these scripts will look for this config file in `~/psp_production.cfg`. So if you want to never think about the config file again, just copy `psp_production.cfg`, which is in the top-directory of the PSP repo, to your home directory (i.e. something like `/Users/lev`). Otherwise, you'll have to supply it explicitly with the --psp_config_path argument.

## Examples

### Use Case 1: dry

You have downloaded an unprocessed gct from Panorama and want to be able to process it (in other words, go from level 2 to level 3 data). You need to use dry.py for this. Navigate to the `broadinstitute_psp` directory of the PSP repository, and type the following in your command line:

```
python dry/dry.py --in_gct_path /Users/lev/Downloads/unprocessed.gct
```

Two output files will be saved to your current directory: the GCT file containing the actual data and a pw file containing various QC things.

### Use Case 2: steep

You have a gct file of new data that you received from a collaborator. You are interested in seeing how it compares to a plate of data that you've made in the past. One way you can do this is by computing the Spearman correlation between each sample of new data with each sample of old data. You need to use steep.py for this. Navigate to the `broadinstitute_psp` directory of the PSP repository, and type the following in your command line:

```
python steep/steep.py --in_gct_path /Users/lev/Downloads/new_data.gct --in_gct2_path /Users/lev/old_data.gct
```

One output file will be saved to your current directory: the GCT file containing the Spearman correlations between the samples of new_data.gct and the samples of old_data.gct.

### Use Case 3: sip

You have the similarity matrix from Use Case 2. Now, you want to know how these similarities compare to all the other stuff that old_data.gct has been compared to before; in other words, you want to compare your similarities to similarities within a reference dataset. You need to use sip.py for this. Navigate to the `broadinstitute_psp` directory of the PSP repository, and type the following in your command line:

```
python sip/sip.py --test_gct_path ./steep_output.gct --bg_gct_path /Users/lev/build_similarity_matrix.gct
```

One output file will be saved to your current directory: the GCT file containing the _connectivities_ between the samples of new_data.gct and the samples of old_data.gct.

Note that an alternative to running steep and sip from the command line is to use the "proteomics-query" web app available on Clue. To query your own P100 or GCP data against Touchstone-P (our reference dataset), visit  [clue.io/proteomics-query](https://clue.io/proteomics-query "Proteomics Query").

Components
----------
- harvest: pushing and pulling data from Panorama (coming soon)
- dry: level 2 &rarr; level 3 data; performs QC
- tear: level 3 &rarr; level 4 data; performs row median normalization or z-scoring    
- steep: computes similarities using level 3 or level 4 data
- sip: computes connectivities using similarity matrices
- external_query: computes connectivities using level 3 or level 4 data (combines steep and sip)
- clue: Python scripts supporting the proteomics-query app on Clue.io
- utils: miscellanous other scripts

Data levels
-----------
![alt text][logo]

[logo]: https://github.com/cmap/psp/blob/a9a5e590757143f09f919e956fbaa206f2be7f52/broadinstitute_psp/misc/2017-07-27_proteomics_data_levels.png "Proteomics Data Levels"
