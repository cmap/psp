#!/bin/bash

source activate psp
cd /usr/src/psp/proteomics-signature-pipeline
python setup_psp.py develop
source deactivate

