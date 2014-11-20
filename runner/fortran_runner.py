import io
import sys, os, errno
import subprocess
import math
import re
import csv
import ConfigParser
import shutil

def read_parameters(args):
    """Parsing  parameters configuration file"""

    with io.open(args.configuration_file, 'r') as configuration_file:
        config = ConfigParser.SafeConfigParser()
        config.readfp(configuration_file)

        parameters = {
            # we only have 2 timesteps from the reference run so we can't validate
            # more timesteps than that
            "NUMBER_OF_TIMESTEPS": 2 if args.validate else config.getint("RUN","NUMBER_OF_TIMESTEPS"),
            "SAVE_INTERVAL": config.getint("RUN","SAVE_INTERVAL"),
            "TIMESTEP": config.getfloat("SIMULATION","TIMESTEP"),
            "SLENDERNESS": config.getfloat("SIMULATION","SLENDERNESS"),
            "NUMBER_OF_TERMS_IN_FORCE_EXPANSION": config.getint("SIMULATION","NUMBER_OF_TERMS_IN_FORCE_EXPANSION"),
            "NUMBER_OF_QUADRATURE_INTERVALS": config.getint("SIMULATION","NUMBER_OF_QUADRATURE_INTERVALS"),
            "NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL": config.getint("SIMULATION","NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL"),
            "GMRES_RESTART": config.getint("GMRES","RESTART"),
            "GMRES_MAX_ITERATIONS": config.getint("GMRES","MAX_ITERATIONS"),
            "GMRES_TOLERANCE": config.getfloat("GMRES","TOLERANCE"),
            "BICGSTAB_MAX_ITERATIONS": config.getint("BICGSTAB","MAX_ITERATIONS"),
            "BICGSTAB_TOLERANCE": config.getfloat("BICGSTAB","TOLERANCE"),
        }

    # Determining number of fibers
    with io.open(args.fibers_file,'r') as fibers:
        parameters['NUMBER_OF_FIBERS'] = int(fibers.readline())

    return parameters

def write_parameters(args, parameters):
    """writes parameters"""

    constants_path = 'fortran/constants.incl'
    with io.open(constants_path, 'w') as constants:
        constants.write(u'#ifndef FIBERS_CONSTANTS_\n')
        constants.write(u'#define FIBERS_CONSTANTS_\n')
        constants.write(u'\n')

        if args.benchmark:
            constants.write(u'#define BENCHMARK\n')
        elif args.validate:
            constants.write(u'#define VALIDATE\n')
        constants.write(u'\n')

        if args.direct:
            constants.write(u'#define DIRECT\n')
        elif args.gmres:
            constants.write(u'#define GMRES\n')
        constants.write(u'\n')

        if args.numerical:
            constants.write(u'#define NUMERICAL\n')
        elif args.analytical:
            constants.write(u'#define ANALYTICAL\n')
        constants.write(u'\n')

        constants.write(u'#define DIMENSIONS (3)\n')
        constants.write(u'#define NUMBER_OF_FIBERS ('+str(parameters['NUMBER_OF_FIBERS'])+')\n')
        constants.write(u'#define TIMESTEP ('+str(parameters['TIMESTEP'])+')\n')
        constants.write(u'#define NUMBER_OF_TIMESTEPS ('+str(parameters['NUMBER_OF_TIMESTEPS'])+')\n')
        constants.write(u'#define SLENDERNESS ('+str(parameters['SLENDERNESS'])+')\n')
        constants.write(u'#define NUMBER_OF_TERMS_IN_FORCE_EXPANSION ('+str(parameters['NUMBER_OF_TERMS_IN_FORCE_EXPANSION'])+')\n')
        constants.write(u'#define NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL ('+str(parameters['NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL'])+')\n')
        constants.write(u'#define NUMBER_OF_QUADRATURE_INTERVALS ('+str(parameters['NUMBER_OF_QUADRATURE_INTERVALS'])+')\n')
        constants.write(u'#define TOTAL_NUMBER_OF_QUADRATURE_POINTS (NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL * NUMBER_OF_QUADRATURE_INTERVALS)\n')
        constants.write(u'\n')

        constants.write(u'#define TOTAL_NUMBER_OF_ROWS (NUMBER_OF_FIBERS * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS)\n')
        constants.write(u'\n')

        constants.write(u'#define GMRES_RESTART ('+str(parameters['GMRES_RESTART'])+')\n')
        constants.write(u'#define GMRES_MAX_ITERATIONS ('+str(parameters['GMRES_MAX_ITERATIONS'])+')\n')
        constants.write(u'#define GMRES_TOLERANCE ('+str(parameters['GMRES_TOLERANCE'])+')\n')

        constants.write(u'#define GMRES_LWORK (GMRES_RESTART * GMRES_RESTART + GMRES_RESTART * (TOTAL_NUMBER_OF_ROWS + 5) + 5 * TOTAL_NUMBER_OF_ROWS + 2)\n')
        constants.write(u'\n')

        constants.write(u'#endif\n')
        constants.write(u'\n')
        constants.flush()
        os.fsync(constants.fileno())

def build(args):
    """Build code"""

    build_path = 'fortran/build/'

    shutil.rmtree(build_path)

    try:
            os.makedirs(build_path)
    except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(build_path):
                    pass
            else:
                    raise

    FNULL = open(os.devnull, 'w')
    if not args.verbose:
        cmake = subprocess.Popen(['cmake', '../'], cwd=build_path, stdout=FNULL)
    else:
        cmake = subprocess.Popen(['cmake', '../'], cwd=build_path)

    if cmake.wait():
        FNULL.close()
        raise Exception("Error building (cmake) the program")

    if not args.verbose:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path, stdout=FNULL)
    else:
        make = subprocess.Popen(['make', '-j8'], cwd=build_path)

    if make.wait():
        FNULL.close()
        raise Exception("Error building (make) the program")

    FNULL.close()

    return build_path

def run(args):
    parameters = read_parameters(args)
    write_parameters(args, parameters)
    build_path = build(args)

    fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers_file])
    if fibers.wait():
        raise Exception("Error running the program")

    print '**************************************************'

def validate(args):
    parameters = read_parameters(args)
    write_parameters(args, parameters)
    build_path = build(args)

    fibers = subprocess.Popen(["./fibers", '../../../'+args.fibers_file], cwd=os.path.join(build_path,'bin/'))
    if fibers.wait():
        raise Exception("Error running the program for validation")

    for i in xrange(parameters['NUMBER_OF_TIMESTEPS']):
        current_a_matrix_filename = os.path.join(build_path,'bin/'+str(i)+'_AMat.out')
        reference_a_matrix_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_AMat.ref'
        current_b_vector_filename = os.path.join(build_path,'bin/'+str(i)+'_BVec.out')
        reference_b_vector_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_BVec.ref'
        current_x_vector_filename = os.path.join(build_path,'bin/'+str(i)+'_XVec.out')
        reference_x_vector_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_XVec.ref'

        current_t_velocity_filename = os.path.join(build_path,'bin/'+str(i)+'_TRANSVel.out')
        reference_t_velocity_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_TRANSVel.ref'
        current_r_velocity_filename = os.path.join(build_path,'bin/'+str(i)+'_ROTVel.out')
        reference_r_velocity_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_ROTVel.ref'

        current_positions_filename = os.path.join(build_path,'bin/'+str(i)+'_POS.out')
        reference_positions_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_POS.ref'
        current_orientations_filename = os.path.join(build_path,'bin/'+str(i)+'_ORIENT.out')
        reference_orientations_filename = 'tests/reference/100_numeric_direct/'+str(i)+'_ORIENT.ref'

        validate = subprocess.Popen(['pypy','tools/validate_mapping.py', 'tests/reference/reference.map', current_a_matrix_filename, 'tests/reference/reference.map', reference_a_matrix_filename])
        if validate.wait():
            raise Exception("Error validating A matrix")
        validate = subprocess.Popen(['pypy','tools/validate_mapping.py', 'tests/reference/reference.map', current_b_vector_filename, 'tests/reference/reference.map', reference_b_vector_filename])
        if validate.wait():
            raise Exception("Error validating B vector")
        validate = subprocess.Popen(['pypy','tools/validate_mapping.py', 'tests/reference/reference.map', current_x_vector_filename, 'tests/reference/reference.map', reference_x_vector_filename])
        if validate.wait():
            raise Exception("Error validating X vector")

        validate = subprocess.Popen(['pypy','tools/validate.py', current_t_velocity_filename, reference_t_velocity_filename])
        if validate.wait():
            raise Exception("Error validating translational velocity")
        validate = subprocess.Popen(['pypy','tools/validate.py', current_r_velocity_filename, reference_r_velocity_filename])
        if validate.wait():
            raise Exception("Error validating rotational velocity")

        validate = subprocess.Popen(['pypy','tools/validate.py', current_positions_filename, reference_positions_filename])
        if validate.wait():
            raise Exception("Error validating positions")
        validate = subprocess.Popen(['pypy','tools/validate.py', current_orientations_filename, reference_orientations_filename])
        if validate.wait():
            raise Exception("Error validating orientations")

def benchmark(args):
    FNULL = open(os.devnull, 'w')
    results = {}

    tests= []
    for i in xrange(1,2048+1):
        if i % 32 == 0:
            tests.append(i)
        elif i % 100 == 0:
            tests.append(i)

    tests = tests[:30]
    tests.reverse()

    for idx,number_of_fibers in enumerate(tests):

        print '  [BENCHMARK]   : ' + str(number_of_fibers) + ' fibers ('+str(idx+1)+'/'+str(len(tests))+')'

        iterations = 4
        benchmark = []

        # allocate dict holding the row of data belonging to the current
        # number of fibers
        results[number_of_fibers] = {}

        sample_mean = 0.0
        sample_deviation = 0.0
        standard_error = 0.0
        relative_standard_error = sys.float_info.max

        while relative_standard_error > args.max_rse:
            for i in xrange(iterations):
                print '                : iteration ('+str(i+1)+'/'+str(iterations)+')'

                gen = subprocess.Popen(['./tools/gen', str(number_of_fibers)], stdout=FNULL)
                if gen.wait():
                    raise Exception("Error generating fibers")
                scene = 'XcT_gen'+str(number_of_fibers)+'.in'

                args.fibers_file = scene
                parameters = read_parameters(args)
                write_parameters(args, parameters)
                build_path = build(args)

                fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), scene], stdout=subprocess.PIPE)

                total = 0.0
                count = 0

                data = {}

                for line in fibers.stdout:
                    line = line.strip()
                    if line.startswith('BENCHMARK'):
                        line = [x.strip() for x in line.split(':')]
                        step_name = str(line[1].strip()).upper()
                        step_value = float(line[2])

                        if not step_name in data:
                            data[step_name] = []
                        data[step_name].append(step_value)

                for step in data:
                    # ignore first timing as warmup
                    times = data[step][1:]
                    data[step] = sum(times)/len(times)

                benchmark.append(data)
                os.remove(scene)

                if fibers.wait():
                    FNULL.close()
                    raise Exception("Error running the program for benchmarking")

            # reset the mean value for the different steps
            for run in benchmark:
                for step in run.keys():
                    results[number_of_fibers][step] = 0.0

            run_sum = reduce(lambda memo, run: memo + run['TOTAL'], benchmark, 0.0)

            sample_mean = run_sum/len(benchmark)
            sample_deviation = 0.0
            for idx, run in enumerate(benchmark):
                sample_deviation += (run['TOTAL'] - sample_mean)**2

                # calculate cumulative moving average
                for step in run.keys():
                    results[number_of_fibers][step] = results[number_of_fibers][step] + (run[step] - results[number_of_fibers][step]) / (idx+1)

            sample_deviation /= len(benchmark)-1
            sample_deviation = math.sqrt(sample_deviation)

            standard_error = sample_deviation / math.sqrt(len(benchmark))
            relative_standard_error = standard_error / sample_mean

            if relative_standard_error > args.max_rse:
                iterations = len(benchmark)
                print '                : Relative Standard Error: ' + str(round(relative_standard_error*100)) + '% - increasing iterations to ' + str(iterations)

        results[number_of_fibers]['TOTAL'] = sample_mean
        results[number_of_fibers]['TOTAL_STD'] = sample_deviation

    FNULL.close()

    with open("results.csv", "wb") as csvfile:
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

    print '**************************************************'
