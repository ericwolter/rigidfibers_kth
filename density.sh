#!/bin/bash

cd tools/
./gen sphere 2000 655.36
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 327.68
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 163.84
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 81.92
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 40.96
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 20.48
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 10.24
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 5.12
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 2.56
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 1.28
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.64
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.32
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.16
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.08
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.04
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.02
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in

cd tools/
./gen sphere 2000 0.01
cd ..
./fibers_runner.py fortran run --gmres --numerical tests/config.ini tools/XcT_gen2000.in
