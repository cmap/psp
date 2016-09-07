# Proteomics Signature Pipeline (PSP)

This is a package of Python scripts that enable processing and analysis of proteomic signatures produced by the LINCS Proteomic Characterization Center for Signaling and Epigenetics (PCCSE) at the Broad Institute. You can download the raw data from the [Panorama Repository](https://panoramaweb.org/labkey/project/LINCS/begin.view? "Panorama Repository"). (You will want the unprocessed gcts.)  
  
A dependency of PSP is the [l1ktools repository](https://github.com/cmap/l1ktools "l1ktools"), a public repo maintained by the Connectivity Map at the Broad Institute. This is where we get the files for parsing, writing, and concatenating gct and gctx files.

## Maintainer

Lev Litichevskiy  
lev@broadinstitute.org  
August 2016

## Setting up your environment

### To set up your environment the first time:

  1. First, you must clone this repository into a local directory on your computer. For example, I cloned the repository in the directory `/Users/lev/code/proteomics-signature-pipeline`. If you need more information about cloning, go to this page provided by Github: https://help.github.com/articles/cloning-a-repository/.

  2. To manage our Python environment, we'll use a program called conda. Download conda from the following website: http://conda.pydata.org/miniconda.html. Miniconda or Anaconda will do, but I'd recommend Miniconda because it's more lightweight.

  3. Now, we will continue with the setup using the command line, so open the terminal or command prompt.

  4. Type `conda info` to verify that conda has been installed on your system. You should see some information about the "Current conda install." If you don't, then conda has not been installed properly.

  5. We will now create an environment with conda that will allow us to use PSP. Type the following in your Terminal:

      ```
      conda create --name psp_env python=2 numpy scipy pandas pytables
      ```
      
      'psp_env' will be the name of your conda environment, and the things after it are the packages that you want that environment to contain. Note that we are using python2, rather than python3. Click 'yes' through the various installation steps.

  6. To activate your environment, type `source activate psp_env`, or if you are on a Windows computer, `activate psp_env`. You should now see `[psp_env]` or `(psp_env)` prepended to the start of your command prompt. For example:

      ```
      (psp_env) /Users/lev/code/proteomics-signature-pipeline $
      ```

  7. We will now run a script to make our environment aware of the contents of the PSP repository that we cloned. Make sure you are in the directory where you cloned your repository (e.g. `/Users/lev/code/proteomics-signature-pipeline`), and then type:

      ```
      python setup_psp.py develop
      ```

  8. This repository is dependent on another repository ([l1ktools](https://github.com/cmap/l1ktools "l1ktools")) for the code that parses and writes gct and gctx files. So you must also clone l1ktools onto your computer (like step 1). For example, I cloned l1ktools into `/Users/lev/code/l1ktools`.   
  
  9. We will now a run a script to make our environment aware of files that we need. (You may need to checkout a different branch first: try `git checkout add_GCToo`.) Navigate to the `io` folder and run another setup script by typing the following:
    
      ```
      cd python/cmap/io
      python setup_GCToo.py develop
      ```
  
  10. That's it! To make sure that everything has been set up correctly, navigate back to top directory of PSP and try executing one of the Python test scripts:

      ```
      cd /Users/lev/code/proteomics-signature-pipeline
      python dry/test_dry.py
      ```
      
      Note that PSP can be run from anywhere. So you could run Python scripts by supplying their full paths. For example:
      
      ```
      python /Users/lev/code/proteomics-signature-pipeline/dry/test_dry.py
      ```
      
      However, before doing this, make sure to check out the section below about the configuration file.
  
### To set up your environment after the first time:

  1. Activate your conda environment by typing `source activate psp_env`, or if you are on a Windows computer, `activate psp_env`.
  2. It's easiest to run your code from the top directory of the PSP repo, so navigate to that directory (the one that contains the dry, steep, etc. directories). For example,
  
    ```
    cd /Users/lev/code/proteomics-signature-pipeline
    ```
    
  3. To make sure again that everything is set up correctly, try executing one of the test scripts again:

      ```
      python dry/test_dry.py
      ```

## Configuration file

Several of the files in PSP -- such as dry, steep, and sip -- require a configuration file. This file contains some parameters used by the production code (e.g. what values to consider NaN, default thresholds, provenance code abbreviations, etc.). By default, these scripts will look for this config file in `~/psp_production.cfg`. So if you want to never think about the config file again, just copy `psp_production.cfg`, which is in the top-directory of the PSP repo, to your home directory. Otherwise, you'll have to supply it explicitly with the --psp_config_path argument.

## Examples

### Use Case 1: dry

You have downloaded an unprocessed gct from Panorama and want to be able to process it (in other words, go from level 2 to level 3 data). You need to use dry.py for this. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python dry/dry.py --in_gct_path /Users/lev/Downloads/unprocessed.gct
```

Two output files will be saved to your current directory: the GCT file containing the actual data and a pw file containing various QC things.

### Use Case 2: steep

You have a gct file of new data that you received from a collaborator. You are interested in seeing how it compares to a plate of data that you've made in the past. One way you can do this is by computing the Spearman correlation between each sample of new data with each sample of old data. You need to use steep.py for this. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python steep/steep.py --in_gct_path /Users/lev/Downloads/new_data.gct --in_gct2_path /Users/lev/old_data.gct
```

One output file will be saved to your current directory: the GCT file containing the Spearman correlations between the samples of new_data.gct and the samples of old_data.gct.

### Use Case 3: sip

You have the similarity matrix from Use Case 2. Now, you want to know how these similarities compare to all the other stuff that old_data.gct has been compared to before; in other words, you want to compare your similarity to the similarity matrix of the "build" or the "corpus." You need to use sip.py for this. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python sip/sip.py --test_gct_path ./steep_output.gct --bg_gct_path /Users/lev/build_similarity_matrix.gct
```

One output file will be saved to your current directory: the GCT file containing the _connectivities_ between the samples of new_data.gct and the samples of old_data.gct.

### Use Case 4: concat_gctoo

A. You have a bunch of files that start with 'LINCS_GCP' in your Downloads folder that you want to concatenate. Type the following in your command line:

```
python /Users/lev/code/l1ktools/python/cmap/io/GCToo/concat_gctoo.py --file_wildcard '/Users/lev/Downloads/LINCS_GCP*'
```

This will save a file called `concated.gct` in your current directory.  

B. You have 2 files that you want to concatenate: /Users/lev/file_to_concatenate1.gct and /Users/lev/file_to_concatenate2.gct. Type the following in your command line:

```
python /Users/lev/code/l1ktools/python/cmap/io/GCToo/concat_gctoo.py --list_of_gct_paths /Users/lev/file_to_concatenate1.gct /Users/lev/file_to_concatenate2.gct
```

C. You have 2 GCToo objects in memory that you want to concatenate. hstack is the method in concat_gctoo.py that actually does the concatenation. From within the Python console or script where you have your 2 GCToos (gct1 & gct2), type the following:

```
import cmap.io.GCToo.concat_gctoo as cg
concated = cg.hstack([gct1, gct2])
```

Components
----------
harvest: pushing and pulling data from Panorama (coming soon)  
dry: level 2 &rarr; level 3 data; performs filtering and normalization  
steep: level 3 or level 4 (z-scored) &rarr; level 5 data; computes similarity  
sip: level 5 &rarr; level 6 data; computes connectivity  
utils: miscellanous other scripts  

Data levels
-----------
![alt text][logo]

[logo]: https://github.com/cmap/proteomics-signature-pipeline/blob/1907ca5661ae617e03678e2e800f06b5503b4b29/2016-07-29_proteomics_data_levels.png "Proteomics Data Levels"
