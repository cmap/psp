# Proteomics Signature Pipeline (PSP)

This is a package of Python scripts that enable processing and analysis of proteomic
signatures produced by the LINCS Proteomic Characterization Center
for Signaling and Epigenetics (PCCSE) at the Broad Institute. You
can download the raw data from the [Panorama Repository](https://panoramaweb.org/labkey/project/LINCS/begin.view? "Panorama Repository"). (You will want the unprocessed gcts.)  

## Maintainer

Lev Litichevskiy  
lev@broadinstitute.org  
August 2016

## Setting up your environment

  1. First, you must clone this repository into a local directory on your computer. For example, I cloned the repository in the directory `/Users/lev/code/proteomics-signature-pipeline`. If you need more information about cloning, go to this page provided by Github: https://help.github.com/articles/cloning-a-repository/.

  2. To manage our Python environment, we'll use a program called conda. Download conda from the following website: http://conda.pydata.org/miniconda.html. Miniconda or Anaconda will do, but I'd recommend Miniconda because it's more lightweight.

  3. Now, we will continue with the setup using the command line, so open the terminal or command prompt.

  4. Type `conda info` to verify that conda has been installed on your system. You should see some information about the "Current conda install." If you don't, then conda has not been installed properly.

  5. We will now create an environment with conda that will allow us to use PSP. Type the following in your Terminal:

  ```
  conda create --name psp_env python=2 numpy scipy pandas pytables
  ```
  
  'psp_env' will be the name of your conda environment, and the things after it are the packages that you want that environment
  to contain. Note that we are using python2, rather than python3. Click 'yes' through the various installation steps.

  6. To activate your environment, type `source activate psp_env`, or if you are on a Windows computer, `activate psp_env`. You should
now see `[psp_env]` or `(psp_env)` prepended to the start of your command prompt. For example:

  ```
  (psp_env) /Users/lev/code/proteomics-signature-pipeline $
  ```

  7. Finally, we will run one more script to make our environment aware of the contents of the PSP repository that we cloned. Make sure you are in the directory where you cloned your repository, and then type:

  ```
  python setup.py develop
  ```

  8. To make sure that everything has been set up correctly, try executing one of the Python test scripts:

  ```
  python dry/test_dry.py
  ```

## Examples

### Use Case 1: dry

You have downloaded an unprocessed gct from Panorama and want to be able to process it (in other words, go from level 2 to level 3 data). You need to use dry.py for this. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python dry/dry.py unprocessed.gct /path/to/output/dir/
```

Two output files will be saved to /path/to/output/dir/: the GCT file containing the actual data and a pw file containing various QC things.

### Use Case 2: steep

You have a gct file of new data that you received from a collaborator. You are interested in seeing how it compares to a plate of data that you've made in the past. One way you can do this is by computing the Spearman correlation between each sample of new data with each sample of old data. You need to use steep.py for this. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python steep/steep.py new.gct /path/to/output/dir/ new_v_old.gct -in_gct2 old.gct
```

One output file will be saved to /path/to/output/dir/: the GCT file containing the Spearman correlations between the samples of new.gct and the samples of old.gct.

### Use Case 3: concat_gctoo

  A) You have a bunch of files that start with 'LINCS_GCP' in /some/dir that you want to concatenate. You want to save the output gct as concated.gct. Navigate to the top directory of the PSP repository, and type the following in your command line:

```
python in_out/concat_gctoo.py '/some/dir/LINCS_GCP*' concated.gct
```

  B) You have 2 GCToo objects in memory that you want to concatenate. hstack is the method in concat_gctoo.py that actually does the concatenation. From within the Python console or script where you have your 2 GCToos (gct1 & gct2), type the following:

```
import in_out.concat_gctoo as cg
concated = cg.hstack([gct1, gct2])
```

Components
----------
harvest: coming soon  
dry: level 2 &rarr; level 3 data; performs filtering and normalization  
steep: computes similarities  
in_out: input/output & concatenation scripts  
utils: miscellanous other scripts  

Data levels
-----------
![alt text][logo]

[logo]: https://github.com/cmap/proteomics-signature-pipeline/blob/1907ca5661ae617e03678e2e800f06b5503b4b29/2016-07-29_proteomics_data_levels.png "Proteomics Data Levels"

More info
---------
Please see the [spec for PSP](https://docs.google.com/a/broadinstitute.com/document/d/1A6-q4ss4JuP-pDkBKMpnCvA2C4KT6JaSxlv6eX2fnx4/edit?usp=sharing "Spec for PSP").
