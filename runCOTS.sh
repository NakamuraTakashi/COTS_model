#!/bin/sh
#
rm *.exe
# Runscript for COTS model:
#
# ------------------ COMPILE & RUN ----------------------------
#
#gfortran -fbounds-check  mod_netcdf.F90 mod_p_coral.F90 mod_cots.F90 main_cots.F90 -I/usr/include -L/usr/lib -lnetcdff -O2 -o cots.exe

gfortran -fbounds-check mod_netcdf.F90 mod_p_coral.F90 mod_cots.F90 mod_cots_larvae.F90 main_cots2.F90 -I/usr/include -L/usr/lib -lnetcdff -fopenmp -O2 -o cots.exe
rm *.mod

export OMP_NUM_THREADS=12
#
./cots.exe
#
