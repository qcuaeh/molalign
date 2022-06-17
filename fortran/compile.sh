#!/bin/sh
gfortran -fPIC -c messages.f90
gfortran -fPIC -c lapack.f90
gfortran -fPIC -c sorting.f90
gfortran -fPIC -c hungarian.f90
gfortran -fPIC -c utilities.f90
gfortran -fPIC -c printing.f90
gfortran -fPIC -c chemistry.f90
gfortran -fPIC -c globals.f90
gfortran -fPIC -c random.f90
gfortran -fPIC -c decoding.f90
gfortran -fPIC -c assignment.f90
gfortran -fPIC -c assortment.f90
gfortran -fPIC -c translation.f90
gfortran -fPIC -c rotation.f90
gfortran -fPIC -c alignment.f90
gfortran -fPIC -c readwrite.f90
gfortran -fPIC -c biasing-r.f90
gfortran -fPIC -c remapping.f90
#gfortran -fPIC -c main.f90
f2py3.4 -c -m ralign ralign.f90 *.o -llapack -L/usr/lib64/atlas
