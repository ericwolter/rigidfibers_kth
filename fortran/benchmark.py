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

tests = tests[:30]
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

    sample_mean = 0.0
    sample_deviation = 0.0
    standard_error = 0.0
    relative_standard_error = sys.float_info.max

    while relative_standard_error > 100000.0:
        for i in xrange(iterations):
            run = open('run100_iter.in', 'r')
            fibers = subprocess.Popen('./ADVECT_FIBERS', stdin=run, stdout=subprocess.PIPE)

            total = 0.0
            count = 0
            for line in fibers.stdout:
                line = line.strip()
                if line.startswith('BENCHMARK'):
                    count += 1
                    if count > 1: # ignore first step as warmup
                        total += float(line.split(':')[1])
            benchmark.append(total/(count-1))

            if fibers.wait():
                run.close()
                FNULL.close()
                raise Exception("Error running the program for benchmarking")

            run.close()

        sample_mean = sum(benchmark)/len(benchmark)
        sample_deviation = 0.0
        for x in benchmark:
            sample_deviation += (x - sample_mean)**2
        sample_deviation /= len(benchmark)-1
        sample_deviation = math.sqrt(sample_deviation)

        standard_error = sample_deviation / math.sqrt(len(benchmark))
        relative_standard_error = standard_error / sample_mean

        iterations = len(benchmark)

    results[number_of_fibers] = sample_mean

FNULL.close()

with open("results.csv", "wb") as csvfile:
    resultswriter = csv.writer(csvfile, dialect="excel-tab")
    row = ["X","SINGLE_AVG_TIME","MULTI_AVG_TIME"]
    resultswriter.writerow(row)
    print '**************************************************'
    print 'Benchmark:'
    print '  '+'\t'.join([str(x) for x in row])
    for key in sorted(results):
        row = [key,results[key],results[key]/8.0]
        resultswriter.writerow(row)
        print '  '+'\t'.join([str(x) for x in row])

print '**************************************************'
