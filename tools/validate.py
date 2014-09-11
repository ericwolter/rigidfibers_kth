import sys
import math

f_path = sys.argv[1]
c_path = sys.argv[2]

print f_path
print c_path

f = [line.strip().split() for line in open(f_path)]
c = [line.strip().split() for line in open(c_path)]

count = 0.0
sum = 0.0
max = -1.0
max_row = 0
max_col = 0
for row in xrange(len(f)):
    f_row = f[row]
    c_row = c[row]
    for col in xrange(len(f_row)):
        count+= 1
        diff = abs(float(f_row[col]) - float(c_row[col]))
        sum += diff
        if diff > max:
            max = diff
            max_row = row
            max_col = col

print 'count', count
print 'sum', sum
print 'max', max, '(',max_row,',',max_col,')'
print 'max_f', f[max_row][max_col]
print 'max_c', c[max_row][max_col]
print 'avg', sum / count