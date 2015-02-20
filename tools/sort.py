from random import shuffle

fibers = []

with open('XcT_gen2000.in') as f:
    lines = [line.rstrip('\n') for line in f]

    M = int(lines[0])
    lines = lines[1:]

    for m in xrange(M):
        p = [float(x) for x in lines[m*2+0].split()]
        o = [float(x) for x in lines[m*2+1].split()]
        fibers.append((p,o))

#random
#shuffle(fibers)
#export = open('XcT_gen'+str(M)+'_random.in', 'w')

#zorder
fibers.sort(key=lambda f: f[0][2],reverse=True)
export = open('XcT_gen'+str(M)+'_zorder.in', 'w')

export.write(str(M) + '\n')
for i in xrange(M):
    f = fibers[i]
    p = f[0]
    o = f[1]
    export.write('\t'.join([str(x) for x in p]) + '\n')
    export.write('\t'.join([str(x) for x in o]) + '\n')

export.close()
