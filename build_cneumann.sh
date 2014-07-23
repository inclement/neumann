#!/bin/bash

# Builds the cyknot module

# Tell it where the numpy header files are
export CFLAGS=-I/usr/lib/python2.7/site-packages/numpy/core/include/ 

python2 setup.py build_ext --inplace
