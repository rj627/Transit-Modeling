#!/bin/bash
echo "This script runs extractions for all available transit observations for WASP-52, WASP-62, WASP-63, WASP-67, WASP-79, WASP-101."
echo "Beginning all extractions. Working on extracting WASP-52, Channel 1"
python spitzer_extraction.py WASP_52 'WASP-52 b' 1 61517056 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-52 Channel 1. Moving on to extracting WASP-62, Channel 1"
python spitzer_extraction.py WASP_62 'WASP-62 b' 1 62167040 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-62 Channel 1. Moving on to extracting WASP-62, Channel 2"
python spitzer_extraction.py WASP_62 'WASP-62 b' 2 62167552 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-62 Channel 2. Moving on to extracting WASP-63, Channel 1"
python spitzer_extraction.py WASP_63 'WASP-63 b' 1 62168064 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-63 Channel 1. Moving on to extracting WASP-63, Channel 2"
python spitzer_extraction.py WASP_63 'WASP-63 b' 2 62168576 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-63 Channel 2. Moving on to extracting WASP-67, Channel 1"
python spitzer_extraction.py WASP_67 'WASP-67 b' 1 62169088 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-67 Channel 1. Moving on to extracting WASP-67, Channel 2"
python spitzer_extraction.py WASP_67 'WASP-67 b' 2 62169600 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-67 Channel 2. Moving on to extracting WASP-79, Channel 1"
python spitzer_extraction.py WASP_79 'WASP-79 b' 1 62173184 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-79 Channel 1. Moving on to extracting WASP-79, Channel 2"
python spitzer_extraction.py WASP_79 'WASP-79 b' 2 62173696 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-79 Channel 2. Moving on to extracting WASP-101, Channel 1"
python spitzer_extraction.py WASP_101 'WASP-101 b' 1 62158336 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with WASP-101 Channel 1. Moving on to extracting WASP-101, Channel 2"
python spitzer_extraction.py WASP_101 'WASP-101 b' 2 62159360 rahulj.4999@gmail.com rahulj.4999@gmail.com

echo "Done with all! Examine xyb plots for consistency, and then run fitting script."
