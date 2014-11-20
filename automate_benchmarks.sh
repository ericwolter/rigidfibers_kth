#!/bin/bash

# stop on error
# set -e

# echo command
set -x

# FORTRAN
# ====

# Direct Solver
# -------------
# ./fibers_runner.py fortran benchmark --direct --numerical tests/config.ini
# cp -f results.csv docs/paper/benchmarks/openmp_direct_numerical.csv
#
# ./fibers_runner.py fortran benchmark --direct --analytical tests/config.ini
# cp -f results.csv docs/paper/benchmarks/openmp_direct_analytical.csv
#
# # GMRES
# # -------------
# ./fibers_runner.py fortran benchmark --gmres --numerical tests/config.ini
# cp -f results.csv docs/paper/benchmarks/openmp_gmres_numerical.csv
#
# ./fibers_runner.py fortran benchmark --gmres --analytical tests/config.ini
# cp -f results.csv docs/paper/benchmarks/openmp_gmres_analytical.csv
#
# # CUDA
# # ====
#
# # MAGMA
# # -------------
# ./fibers_runner.py cuda benchmark --magma --numerical --D2 tests/config.ini
# cp -f results.csv docs/paper/benchmarks/cuda_magma_numerical_2D.csv
#
# # crashes for unknown reasons...
# # ./fibers_runner.py cuda benchmark --magma --analytical --D2 tests/config.ini
# # cp -f results.csv docs/paper/benchmarks/cuda_magma_analytical_2D.csv
#
# ./fibers_runner.py cuda benchmark --magma --numerical --D1 tests/config.ini
# cp -f results.csv docs/paper/benchmarks/cuda_magma_numerical_1D.csv
#
# ./fibers_runner.py cuda benchmark --magma --numerical --D3 tests/config.ini
# cp -f results.csv docs/paper/benchmarks/cuda_magma_numerical_3D.csv
#
# GMRES
# -------------
./fibers_runner.py cuda benchmark --gmres --numerical --D2 tests/config.ini
cp -f results.csv docs/paper/benchmarks/cuda_gmres_numerical_2D.csv

# BICGSTAB
# -------------
./fibers_runner.py cuda benchmark --bicgstab --numerical --D2 tests/config.ini
cp -f results.csv docs/paper/benchmarks/cuda_bicgstab_numerical_2D.csv
