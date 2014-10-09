#!/usr/bin/python

#####
#
# Parsing command line arguments
#
#####
import argparse
parser = argparse.ArgumentParser(description='Run helper for rigid fibers')
parser.add_argument('config',metavar='CONFIG_FILE',type=open,help='the configuration file')
parser.add_argument('fibers',metavar='FIBER_FILE',type=open,help='the fiber file containing the initial positions and orientations')
parser.add_argument('--benchmark',action='store_true',help='an option to dump configuration file (default: False)')
parser.add_argument('--validate',action='store_true',help='an option to validate results against reference implementation (default: False) (note: limits number of steps to 1)')
parser.add_argument('--force1D',action='store_true',help='an option to force to execute kernels in 1D (default: False)')

args = parser.parse_args()

#####
#
# Parsing  parameters configuration file
#
#####
from decimal import Decimal

parameters = {}
for idx,line in enumerate(args.config):
	line = line.strip()
	if not line:
		continue
	value = line.split()[0]
	if idx==0:
		parameters['label'] = int(value)
	elif idx==1:
		pass
	elif idx==2:
		parameters['number_of_terms_in_force_expansion'] = int(value)
	elif idx==3:
		parameters['slenderness'] = Decimal(value)
	elif idx==4:
		parameters['timestep'] = Decimal(value)
	elif idx==5:
		parameters['number_of_timesteps'] = 1 if args.validate else int(value)
	elif idx==6:
		pass
	elif idx==7:
		pass
	elif idx==8:
		parameters['number_of_quadrature_intervals'] = int(value)
	elif idx==9:
		pass
	elif idx==10:
		pass

args.config.close()

#####
#
# Determining number of fibers
#
#####

parameters['number_of_fibers'] = int(args.fibers.readline())
args.fibers.close()

#####
#
# Write cuda constants kernel
#
#####
import io

cuda_constants_path = 'src/kernels/constants.cu'
with io.open(cuda_constants_path, 'w') as cuda:
	cuda.write(u'#ifndef FIBERS_CONSTANTS_KERNEL_\n')
	cuda.write(u'#define FIBERS_CONSTANTS_KERNEL_\n')
	cuda.write(u'\n')
	
	cuda.write(u'#include "../common.h"\n')
	cuda.write(u'\n')
	
        if args.benchmark:
            cuda.write(u'#define BENCHMARK\n')
        if args.validate:
            cuda.write(u'#define VALIDATE\n')
        if args.force1D:
            cuda.write(u'#define FORCE_1D\n')
        cuda.write(u'\n')

	cuda.write(u'#define DIMENSIONS (3)\n')
	cuda.write(u'#define NUMBER_OF_FIBERS ('+str(parameters['number_of_fibers'])+')\n')
	cuda.write(u'#define TIMESTEP ('+str(parameters['timestep'])+')\n')
	cuda.write(u'#define NUMBER_OF_TIMESTEPS ('+str(parameters['number_of_timesteps'])+')\n')
	cuda.write(u'#define SLENDERNESS ('+str(parameters['slenderness'])+')\n')
	cuda.write(u'#define NUMBER_OF_TERMS_IN_FORCE_EXPANSION ('+str(parameters['number_of_terms_in_force_expansion'])+')\n')
	cuda.write(u'#define NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL (3)\n')
	cuda.write(u'#define NUMBER_OF_QUADRATURE_INTERVALS ('+str(parameters['number_of_quadrature_intervals'])+')\n')
	cuda.write(u'#define TOTAL_NUMBER_OF_QUADRATURE_POINTS (NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL * NUMBER_OF_QUADRATURE_INTERVALS)\n')
        #cuda.write(u'#define USE_ANALYTICAL_INTEGRATION (false)\n')
	cuda.write(u'\n')

	cuda.write(u'#define TOTAL_NUMBER_OF_ROWS (NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS)\n')
	cuda.write(u'\n')

	cuda.write(u'__constant__ float quadrature_points[TOTAL_NUMBER_OF_QUADRATURE_POINTS];\n')
	cuda.write(u'__constant__ float quadrature_weights[TOTAL_NUMBER_OF_QUADRATURE_POINTS];\n')
	cuda.write(u'__constant__ float legendre_polynomials[TOTAL_NUMBER_OF_QUADRATURE_POINTS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION];\n')
	cuda.write(u'__constant__ float lambda[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];\n')
	cuda.write(u'__constant__ float eigen[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];\n')
	cuda.write(u'\n')
	
	cuda.write(u'#endif //FIBERS_CONSTANTS_KERNEL_\n')
	cuda.write(u'\n')
        cuda.flush()

#####
#
# Build code
#
#####
import sys, os, errno
import subprocess

build_path = 'build/'
try:
	os.makedirs(build_path)
except OSError as exc:
	if exc.errno == errno.EEXIST and os.path.isdir(build_path):
		pass
	else:
		raise
	
cmake = subprocess.Popen(['cmake', '../'], cwd=build_path)
if cmake.wait():
	sys.exit(1)

make = subprocess.Popen(['make', '-j8'], cwd=build_path)
if make.wait():
	sys.exit(2)

#####
#
# Run fibers
#
#####
import inspect
fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers.name])
if fibers.wait():
    sys.exit(3)

if args.validate:
    validate = subprocess.Popen(['python','tools/validate.py', os.path.join(build_path,'bin/current.map'), os.path.join(build_path,'bin/a_matrix.out'), os.path.join(build_path,'bin/current.map'), 'tests/reference/100_numeric_gmres/0_AMat.ref'])
    print validate.wait()


