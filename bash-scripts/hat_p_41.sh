#!/bin/bash

echo "Extraction and Fit for HAT_P_41 b"
echo "Extraction"
python spitzer_extraction.py HAT_P_41 'HAT-P-41 b' 1 62152704 rahulj.4999@gmail.com rahulj.4999@gmail.com
python spitzer_extraction.py HAT_P_41 'HAT-P-41 b' 2 62153216 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Fitting"
bash all_fits.sh