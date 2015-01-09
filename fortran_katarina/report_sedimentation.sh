#!/bin/bash

make -j8 && ./ADVECT_FIBERS < run16.in | grep SEDIMENTATION | tr -s ' ' | cut -d ' ' -f 3
