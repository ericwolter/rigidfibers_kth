#!/usr/bin/python

import io
import sys, os, errno
import subprocess
import math
import glob
import re
import csv
from decimal import Decimal

def read_parameters(args):
    """Parsing  parameters configuration file"""

    parameters = {}

    with io.open(args.config,'r') as config:
        for idx,line in enumerate(config):
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

    # Determining number of fibers
    with io.open(args.fibers,'r') as fibers:
        parameters['number_of_fibers'] = int(fibers.readline())

    return parameters

def write_parameters(args, parameters):
    """writes parameters"""

    # Write cuda constants kernel
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
            if args.gmres:
                cuda.write(u'#define GMRES\n')
            if args.magma:
                cuda.write(u'#define MAGMA\n')
            cuda.write(u'\n')

            cuda.write(u'#define DIMENSIONS (3)\n')
            cuda.write(u'#define NUMBER_OF_FIBERS ('+str(parameters['number_of_fibers'])+')\n')
            cuda.write(u'#define TIMESTEP ('+str(parameters['timestep'])+')\n')
            cuda.write(u'#define NUMBER_OF_TIMESTEPS ('+str(parameters['number_of_timesteps'])+')\n')
            cuda.write(u'#define SLENDERNESS ('+str(parameters['slenderness'])+'f)\n')
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

def build(args):
    """Build code"""

    build_path = 'build/'
    try:
            os.makedirs(build_path)
    except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(build_path):
                    pass
            else:
                    raise

    FNULL = open(os.devnull, 'w')
    if args.benchmark:
        cmake = subprocess.Popen(['cmake', '../'], cwd=build_path, stdout=FNULL)
    else:
        cmake = subprocess.Popen(['cmake', '../'], cwd=build_path)

    if cmake.wait():
        FNULL.close()
        sys.exit(1)

    if args.benchmark:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path, stdout=FNULL)
    else:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path)

    if make.wait():
        FNULL.close()
        sys.exit(2)

    FNULL.close()

    return build_path

#####
#
# Parsing command line arguments
#
#####
import argparse
parser = argparse.ArgumentParser(description='Run helper for rigid fibers')
parser.add_argument('config',metavar='CONFIG_FILE',type=str,help='the configuration file')
parser.add_argument('fibers',metavar='FIBER_FILE',type=str,help='the fiber file containing the initial positions and orientations')
parser.add_argument('--benchmark',action='store_true',help='an option to dump configuration file (default: False)')
parser.add_argument('--validate',action='store_true',help='an option to validate results against reference implementation (default: False) (note: limits number of steps to 1)')
parser.add_argument('--force1D',action='store_true',help='an option to force to execute kernels in 1D (default: False)')
parser.add_argument('--gmres',action='store_true',help='an option to use GMRES solver from ViennaCL (default: False)')
parser.add_argument('--magma',action='store_true',help='an option to use direct solver from MAGMA instead of iterative solver from ViennaCL (default: False)')

args = parser.parse_args()

#####
#
# Run fibers
#
#####
if args.validate:
    parameters = read_parameters(args)
    write_parameters(args, parameters)
    build_path = build(args)

    fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers])
    if fibers.wait():
        sys.exit(3)

    validate = subprocess.Popen(['python','tools/validate.py', os.path.join(build_path,'bin/current.map'), os.path.join(build_path,'bin/a_matrix.out'), 'tests/reference/reference.map', 'tests/reference/100_numeric_gmres/0_AMat.ref'])
    validate.wait()
    validate = subprocess.Popen(['python','tools/validate.py', os.path.join(build_path,'bin/current.map'), os.path.join(build_path,'bin/b_vector.out'), 'tests/reference/reference.map', 'tests/reference/100_numeric_gmres/0_BVec.ref'])
    validate.wait()
    validate = subprocess.Popen(['python','tools/validate.py', os.path.join(build_path,'bin/current.map'), os.path.join(build_path,'bin/x_vector.out'), 'tests/reference/reference.map', 'tests/reference/100_numeric_gmres/0_XVec.ref'])
    validate.wait()
elif args.benchmark:

    tests = glob.glob("tests/*.in")

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        """
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        """
        return [atoi(c) for c in re.split('(\d+)', text)]

    tests.sort(key=natural_keys)
    #print tests[:30]

    FNULL = open(os.devnull, 'w')

    results = {}

    for test in tests[:30]:

        print test

        args.fibers = test
        parameters = read_parameters(args)
        write_parameters(args, parameters)
        build_path = build(args)

        iterations = 8
        benchmark = []

        sample_mean = 0.0
        sample_deviation = 0.0
        standard_error = 0.0
        relative_standard_error = sys.float_info.max

        while relative_standard_error > 0.05:
            for i in xrange(iterations):
                fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers], stdout=FNULL)
                if fibers.wait():
                    FNULL.close()
                    sys.exit(3)
                with io.open(os.path.join(build_path,'bin/performance.out'), 'r') as performance:
                    total = 0.0
                    for idx,line in enumerate(performance):
                        if idx > 0: # ignore header
                            line = line.rstrip().split(',')
                            total += float(line[1])
                    benchmark.append(total)

            sample_mean = sum(benchmark)/len(benchmark)
            sample_deviation = 0.0
            for x in benchmark:
                sample_deviation += (x - sample_mean)**2
            sample_deviation /= len(benchmark)-1
            sample_deviation = math.sqrt(sample_deviation)

            standard_error = sample_deviation / math.sqrt(len(benchmark))
            relative_standard_error = standard_error / sample_mean

            #print parameters["number_of_fibers"], iterations, sample_mean, sample_deviation,relative_standard_error
            iterations = len(benchmark)

        results[parameters["number_of_fibers"]] = sample_mean

    print results
    FNULL.close()

    #with open("benchmarks/overall.csv", "wb") as csvfile:
        #resultswriter = csv.writer(csvfile, dialect="excel-tab")
        #resultswriter.writerow(["X","test"])

else:
    parameters = read_parameters(args)
    write_parameters(args, parameters)
    build_path = build(args)

    fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers])
    if fibers.wait():
        sys.exit(3)

print '**************************************************'
