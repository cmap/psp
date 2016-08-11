Proteomics Signature Pipeline (PSP)
=========

This is a package of Python scripts that enable processing and analysis of proteomic
signatures produced by the LINCS Proteomic Characterization Center
for Signaling and Epigenetics (PCCSE) at the Broad Institute. You
can download the raw data from the [Panorama Repository](https://panoramaweb.org/labkey/project/LINCS/begin.view? "Panorama Repository"). (You will want the unprocessed gcts.)  

Maintainer
----------
Lev Litichevskiy  
lev@broadinstitute.org  
August 2016

Setting up your environment
---------------------------
  1. First, you must clone this repository into a local directory on your computer. For example, I cloned the repository in the directory `/Users/lev/code/proteomics-signature-pipeline`. If you need more information about cloning, go to this page provided by Github: https://help.github.com/articles/cloning-a-repository/.

  2. To manage our Python environment, we'll use a program called conda. Download conda from the following website: http://conda.pydata.org/miniconda.html. Miniconda or Anaconda will do, but I'd recommend Miniconda because it's more lightweight.

  3. Now, we will continue with the setup in the Terminal, so open the Terminal.

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

For more information, please refer to the instructions in the [spec for PSP](https://docs.google.com/a/broadinstitute.com/document/d/1A6-q4ss4JuP-pDkBKMpnCvA2C4KT6JaSxlv6eX2fnx4/edit?usp=sharing "Spec for PSP").  

The spec contains several of the most common uses of the PSP repo.

Components
----------
harvest: coming soon  
dry: level 2 -> level 3 data; performs filtering and normalization  
steep: computes similarities
in_out: input/output & concatenation scripts  
utils: miscellanous other scripts  

Data levels
-----------
![alt text][logo]

[logo]: https://github.com/cmap/proteomics-signature-pipeline/blob/1907ca5661ae617e03678e2e800f06b5503b4b29/2016-07-29_proteomics_data_levels.png "Proteomics Data Levels"
