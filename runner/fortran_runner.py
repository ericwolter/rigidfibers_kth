import io
import sys, os, errno
import subprocess
import math
import re
import csv

def build(args):
    """Build code"""

    build_path = 'fortran/build/'
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
    build_path = build(args)
    fibers = subprocess.Popen([os.path.join(build_path,'bin/fibers'), args.fibers_file])
    if fibers.wait():
        raise Exception("Error running the program")

    print '**************************************************'

def validate(args):
    print "FORTRAN validate"
    print args

def benchmark(args):
    FNULL = open(os.devnull, 'w')
    results = {}
    src_path = 'fortran/'

    tests= []
    for i in xrange(1,2048+1):
        if i % 32 == 0:
            tests.append(i)
        elif i % 100 == 0:
            tests.append(i)

    tests = tests[:10]
    tests.reverse()

    for idx,number_of_fibers in enumerate(tests):

        print '  [BENCHMARK]   : ' + str(number_of_fibers) + ' fibers ('+str(idx+1)+'/'+str(len(tests))+')'

        iterations = 4
        benchmark = []

        run = open(os.path.join(src_path,'constants.incl'), 'r')
        rundata = run.read()
        run.close()

        # replace
        rundata = re.sub(r'#define NUMBER_OF_FIBERS \((\d+)\)',"#define NUMBER_OF_FIBERS ("+str(number_of_fibers)+")", rundata)

        run = open(os.path.join(src_path,'constants.incl'), 'w')
        rundata = run.write(rundata)
        run.close()

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

                gen = subprocess.Popen(['python','tools/gen2.py', str(number_of_fibers)], stdout=FNULL)
                if gen.wait():
                    raise Exception("Error generating fibers")
                scene = 'XcT_gen'+str(number_of_fibers)+'.in'

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
                    times = data[step]#[1:]
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

            run_sum = reduce(lambda memo, run: memo + run['$TOTAL'], benchmark, 0.0)

            sample_mean = run_sum/len(benchmark)
            sample_deviation = 0.0
            for idx, run in enumerate(benchmark):
                sample_deviation += (run['$TOTAL'] - sample_mean)**2

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

        results[number_of_fibers]['$TOTAL'] = sample_mean

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
