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
                        # we only have 2 timesteps from the reference run so we can't validate
                        # more timesteps than that
                        parameters['number_of_timesteps'] = 2 if args.validate else int(value)
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
        raise Exception("Error building (cmake) the program")

    if args.benchmark:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path, stdout=FNULL)
    else:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path)

    if make.wait():
        FNULL.close()
        raise Exception("Error building (make) the program")

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
        raise Exception("Error running the program for validation")

    for i in xrange(parameters['number_of_timesteps']):
        current_a_matrix_filename = os.path.join(build_path,'bin/a_matrix_'+str(i)+'.out')
        reference_a_matrix_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_AMat.ref'
        current_b_vector_filename = os.path.join(build_path,'bin/b_vector_'+str(i)+'.out')
        reference_b_vector_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_BVec.ref'
        current_x_vector_filename = os.path.join(build_path,'bin/x_vector_'+str(i)+'.out')
        reference_x_vector_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_XVec.ref'

        current_t_velocity_filename = os.path.join(build_path,'bin/t_vel_'+str(i)+'.out')
        reference_t_velocity_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_TRANSVel.ref'
        current_r_velocity_filename = os.path.join(build_path,'bin/r_vel_'+str(i)+'.out')
        reference_r_velocity_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_ROTVel.ref'

        current_positions_filename = os.path.join(build_path,'bin/positions_'+str(i)+'.out')
        reference_positions_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_POS.ref'
        current_orientations_filename = os.path.join(build_path,'bin/orientations_'+str(i)+'.out')
        reference_orientations_filename = 'tests/reference/100_numeric_gmres/'+str(i)+'_ORIENT.ref'


        validate = subprocess.Popen(['python','tools/validate_mapping.py', os.path.join(build_path,'bin/current.map'), current_a_matrix_filename, 'tests/reference/reference.map', reference_a_matrix_filename])
        if validate.wait():
            raise Exception("Error validating A matrix")
        validate = subprocess.Popen(['python','tools/validate_mapping.py', os.path.join(build_path,'bin/current.map'), current_b_vector_filename, 'tests/reference/reference.map', reference_b_vector_filename])
        if validate.wait():
            raise Exception("Error validating B vector")
        validate = subprocess.Popen(['python','tools/validate_mapping.py', os.path.join(build_path,'bin/current.map'), current_x_vector_filename, 'tests/reference/reference.map', reference_x_vector_filename])
        if validate.wait():
            raise Exception("Error validating X vector")

        validate = subprocess.Popen(['python','tools/validate.py', current_t_velocity_filename, reference_t_velocity_filename])
        if validate.wait():
            raise Exception("Error validating translational velocity")
        validate = subprocess.Popen(['python','tools/validate.py', current_r_velocity_filename, reference_r_velocity_filename])
        if validate.wait():
            raise Exception("Error validating rotational velocity")

        validate = subprocess.Popen(['python','tools/validate.py', current_positions_filename, reference_positions_filename])
        if validate.wait():
            raise Exception("Error validating positions")
        validate = subprocess.Popen(['python','tools/validate.py', current_orientations_filename, reference_orientations_filename])
        if validate.wait():
            raise Exception("Error validating orientations")

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

    def extract_timestep(filename):
        return re.search('_(\d+)\.', filename).groups()[0]

    tests.sort(key=natural_keys)
    #print tests[:30]

    FNULL = open(os.devnull, 'w')

    results = {}

    tests = tests[:30]
    tests.reverse()

    for idx,test in enumerate(tests):

        print '  [BENCHMARK]   : ' + test + ' ('+str(idx+1)+'/'+str(len(tests))+')'

        args.fibers = test
        parameters = read_parameters(args)
        number_of_fibers = parameters["number_of_fibers"]
        write_parameters(args, parameters)
        build_path = build(args)

        iterations = 2
        benchmark = []

        # allocate dict holding the row of data belonging to the current
        # number of fibers
        results[number_of_fibers] = {}

        sample_mean = 0.0
        sample_deviation = 0.0
        standard_error = 0.0
        relative_standard_error = sys.float_info.max

        while relative_standard_error > 0.05:
            for i in xrange(iterations):
                fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers], stdout=FNULL)
                if fibers.wait():
                    FNULL.close()
                    raise Exception("Error running the program for benchmarking")
                performance = os.path.join(build_path,'bin/performance.out')
                with io.open(performance) as performance_file:
                    run = {'$TOTAL':0.0}
                    for idx,line in enumerate(performance_file):
                        if idx > 0: # ignore header
                            line = line.rstrip().split(',')

                            step_name = str(line[0].strip())
                            step_value = float(line[1])

                            run[step_name] = step_value
                            run['$TOTAL'] += step_value

                    benchmark.append(run)
                os.remove(performance)

            # reset the mean value for the different steps
            for run in benchmark:
                for step in run.keys():
                    results[parameters["number_of_fibers"]][step] = 0.0

            run_sum = reduce(lambda memo, run: memo + run['$TOTAL'], benchmark, 0.0)

            sample_mean = run_sum/len(benchmark)
            sample_deviation = 0.0
            for idx, run in enumerate(benchmark):
                sample_deviation += (run['$TOTAL'] - sample_mean)**2

                # calculate cumulative moving average
                for step in run.keys():
                    results[parameters["number_of_fibers"]][step] = results[parameters["number_of_fibers"]][step] + (run[step] - results[parameters["number_of_fibers"]][step]) / (idx+1)

            sample_deviation /= len(benchmark)-1
            sample_deviation = math.sqrt(sample_deviation)

            standard_error = sample_deviation / math.sqrt(len(benchmark))
            relative_standard_error = standard_error / sample_mean

            iterations = len(benchmark)

        results[parameters["number_of_fibers"]]['$TOTAL'] = sample_mean

    FNULL.close()

    with open("benchmarks/results.csv", "wb") as csvfile:
        resultswriter = csv.writer(csvfile, dialect="excel-tab")
        row = ["X"] + [x for x in sorted(results.values()[0].keys())]
        resultswriter.writerow(row)
        print '**************************************************'
        print 'Benchmark:'
        print '  '+'\t'.join([str(x) for x in row])
        for key in sorted(results):
            row = [key] + [results[key][k] for k in sorted(results[key].keys())]
            resultswriter.writerow(row)
            print '  '+'\t'.join([str(x) for x in row])
else:
    parameters = read_parameters(args)
    write_parameters(args, parameters)
    build_path = build(args)

    fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers])
    if fibers.wait():
        raise Exception("Error running the program")

print '**************************************************'
