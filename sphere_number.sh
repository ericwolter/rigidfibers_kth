#!/bin/bash

cd tools/
./gen sphere 2000 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config900.ini tools/XcT_gen2000.in
mkdir /NOBACKUP/ewolter/sync/sphere_2000_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_2000_0040/
rm -f build/bin/*.out

cd tools/
./gen sphere 1500 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config900.ini tools/XcT_gen1500.in
mkdir /NOBACKUP/ewolter/sync/sphere_1500_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_1500_0040/
rm -f build/bin/*.out

cd tools/
./gen sphere 1000 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config600.ini tools/XcT_gen1000.in
mkdir /NOBACKUP/ewolter/sync/sphere_1000_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_1000_0040/
rm -f build/bin/*.out

cd tools/
./gen sphere 500 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config600.ini tools/XcT_gen500.in
mkdir /NOBACKUP/ewolter/sync/sphere_500_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_500_0040/
rm -f build/bin/*.out

cd tools/
./gen sphere 250 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config500.ini tools/XcT_gen250.in
mkdir /NOBACKUP/ewolter/sync/sphere_250_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_250_0040/
rm -f build/bin/*.out

cd tools/
./gen sphere 100 0.4
cd ..
./fibers_runner.py cuda run --magma --numerical --D2 tests/config400.ini tools/XcT_gen100.in
mkdir /NOBACKUP/ewolter/sync/sphere_100_0040/
cp build/bin/*.out /NOBACKUP/ewolter/sync/sphere_100_0040/
rm -f build/bin/*.out


