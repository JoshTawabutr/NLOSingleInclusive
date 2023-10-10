#!/bin/bash
module add CMake
module add GSL/2.6-GCC-9.3.0
g++ main.cpp -L/cvmfs/fgi.csc.fi/apps/el7/GCCcore/7.3.0/bin/gcc -o penmp_oberon -lLHAPDF -O2 -fopenmp `gsl-config --cflags` `gsl-config --libs` `/home/tawabyx/lhapdf/bin/lhapdf-config --cflags` `/home/tawabyx/lhapdf/bin/lhapdf-config --ldflags`
