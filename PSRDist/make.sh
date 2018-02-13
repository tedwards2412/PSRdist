#!/bin/bash

# Script to compile cython version of ymw16

rm -r build/
rm *.so *.c
python setup.py build_ext --inplace
