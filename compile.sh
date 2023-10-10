#!/bin/bash
module add CMake; 
module add GSL/2.5-GCC-7.3.0-2.30;
g++ main.cpp -L/cvmfs/fgi.csc.fi/apps/el7/GCCcore/7.3.0/bin/gcc -o penmp -lLHAPDF -O2 -fopenmp `gsl-config --cflags` `gsl-config --libs` `/home/tawabyx/lhapdf/bin/lhapdf-config --cflags` `/home/tawabyx/lhapdf/bin/lhapdf-config --ldflags`
