Proteomics Signature Pipeline (PSP)
=========

This is a package of Python scripts used for processing proteomic
signatures produced by the LINCS Proteomic Characterization Center
for Signaling and Epigenetics (PCCSE) at the Broad Institute. You
can download the raw data from the [Panorama Repository](https://panoramaweb.org/labkey/project/LINCS/begin.view? "Panorama Repository"). (You will want the unprocessed gcts.)  

![alt text][logo]

[logo]: https://github.com/cmap/proteomics-signature-pipeline/blob/1907ca5661ae617e03678e2e800f06b5503b4b29/2016-07-29_proteomics_data_levels.png "Proteomics Data Levels"

Maintainer
----------
Lev Litichevskiy  
lev@broadinstitute.org  
August 2016

Setting up your environment
---------------------------

The easiest way to set up your environment is by using conda and setup.py.
Please refer to the instructions in the [spec for PSP](https://docs.google.com/a/broadinstitute.com/document/d/1A6-q4ss4JuP-pDkBKMpnCvA2C4KT6JaSxlv6eX2fnx4/edit?usp=sharing "Spec for PSP").

TODO(lev): write the instructions here in add'n to having them in the spec...

Components
----------
harvest: coming soon  
dry: level 2 -> level 3 data; performs filtering and normalization  
steep: level 5 -> level 6 -> level 7 data; computes similarities and connectivities  
in_out: input/output scripts  
utils: miscellanous other scripts  
