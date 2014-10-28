#!/usr/bin/python

import io
import sys, os, errno
import subprocess
import math
import glob
import re
import csv
from decimal import Decimal

def build(args):
    """Build code"""

    build_path = '.'

    FNULL = open(os.devnull, 'w')
    make = subprocess.Popen(['make', '-j8'], cwd=build_path, stdout=FNULL)
    if make.wait():
        FNULL.close()
        raise Exception("Error building (make) the program")
    FNULL.close()

#####
#
# Run fibers
#
#####

build({})

tests = glob.glob("XcT_init*.in")

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    """
    return [atoi(c) for c in re.split('(\d+)', text)]

def extract_timestep(filename):
    return re.search('init(\d+)\.', filename).groups()[0]

tests.sort(key=natural_keys)

FNULL = open(os.devnull, 'w')

results = {}

tests = tests[:10]
tests.reverse()

for idx,test in enumerate(tests):

    print '  [BENCHMARK]   : ' + test + ' ('+str(idx+1)+'/'+str(len(tests))+')'

    iterations = 2
    benchmark = []

    number_of_fibers = int(extract_timestep(test))
    run = open('run100_iter.in', 'r')
    rundata = run.read()
    run.close()

    # replace
    rundata = re.sub(r'(\d+)\s*!Label of indata file XcT_init##',str(number_of_fibers) + " !Label of indata file XcT_init##", rundata)

    run = open('run100_iter.in', 'w')
    rundata = run.write(rundata)
    run.close()

    # allocate dict holding the row of data belonging to the current
    # number of fibers
    results[number_of_fibers] = {}

    sample_mean = 0.0
    sample_deviation = 0.0
    standard_error = 0.0
    relative_standard_error = sys.float_info.max

    while relative_standard_error > 0.05:
        for i in xrange(iterations):
            run = open('run100_iter.in', 'r')
            fibers = subprocess.Popen('./ADVECT_FIBERS', stdin=run, stdout=subprocess.PIPE)


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

            if fibers.wait():
                run.close()
                FNULL.close()
                raise Exception("Error running the program for benchmarking")

            run.close()

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

        iterations = len(benchmark)


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
