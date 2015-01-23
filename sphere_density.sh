#!/bin/bash

#cd tools/
#./gen sphere 2000 2.56
#cd ..
#./fibers_runner.py cuda run --magma --numerical --D2 tests/config2500.ini tools/XcT_gen2000.in
#mkdir /NOBACKUP/ewolter/sync/sphere_2000_0256/
#cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0256/
#rm -f build/bin/*.out

#cd tools/
#./gen sphere 2000 1.92
#cd ..
#./fibers_runner.py cuda run --magma --numerical --D2 tests/config2000.ini tools/XcT_gen2000.in
#mkdir /NOBACKUP/ewolter/sync/sphere_2000_0192/
#cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0192/
#rm -f build/bin/*.out

cd tools/
./gen sphere 2000 1.28
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config1400.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0128/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0128/
rm -f build/bin/*.out

cd tools/
./gen sphere 2000 0.96
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config1200.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0096/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0096/
rm -f build/bin/*.out

cd tools/
./gen sphere 2000 0.64
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config900.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0064/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0064/
rm -f build/bin/*.out

cd tools/
./gen sphere 2000 0.32
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config600.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0032/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0032/
rm -f build/bin/*.out

cd tools/
./gen sphere 2000 0.16
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config500.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0016/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0016/
rm -f build/bin/*.out

cd tools/
./gen sphere 2000 0.08
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config400.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0008/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0008/
rm -f build/bin/*.out


