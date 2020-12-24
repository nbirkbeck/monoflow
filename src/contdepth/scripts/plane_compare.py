#!/bin/env python
# 
# Quick hack to draw some 3D diagrams for the results after
# trying several of the options.
#
import os, sys
import re
from numpy import * 
import mpl_toolkits.mplot3d.axes3d as p3
import pylab as p

class Stats:
    def __init__(self, fimg, vdisp):
        self.fimg = fimg
        self.vdisp = vdisp
        self.seqs = []

    def average(self):
        self.average = [0,0,0]
        for i in xrange(0, len(self.seqs)):
            self.average[0] += self.seqs[i][0]
            self.average[1] += self.seqs[i][1][0]
            self.average[2] += self.seqs[i][2]

        self.average = [a/len(self.seqs) for a in self.average]

def parse_block(lines):
    p = re.match(r'-----fimg:(?P<fimg>[0-9\.]*) vdisp:(?P<vdisp>[0-9\.]*)', 
                 lines[0]);
    fimg = p.groupdict()['fimg']
    vdisp = p.groupdict()['vdisp']
    
    s = Stats(float(fimg), float(vdisp))

    seq = 0
    seqs = []

    for li in xrange(1, len(lines)):
        l = lines[li]
        l = l.strip()

        if l.startswith('sequence:'):
            seq = int(l[9:])
            seqs.append([])

        if l.startswith('disp:'):
            p = eval(l[5:])
            seqs[seq].append(p)

        if l.startswith('image'):
            seqs[seq].append(float(l.split(':')[1].strip().split(' ')[0]))
        
        if l.startswith('----'):
            break

    s.seqs = seqs
    s.average()

    print fimg, vdisp
    print seqs
    print 
    return s, li

def plot(stats, index, zlabel = 'Z'):
    fimg, vdisp = {}, {}

    for s in stats:
        fimg[s.fimg] = 1
        vdisp[s.vdisp] = 1
    
    fimg = fimg.keys()
    vdisp = vdisp.keys()

    fimg.sort()
    vdisp.sort()

    fimg_index = dict(zip(fimg, range(0, len(fimg))))
    vdisp_index = dict(zip(vdisp, range(0, len(vdisp))))
    x = array([0.0]*len(fimg))
    y = array([0.0]*len(vdisp))
    data = outer(x, y)
    x = outer(array(fimg), ones([1, len(vdisp)]))
    y = outer(array(vdisp), ones([1, len(fimg)])).transpose()
    
    for s in stats:
        data[fimg_index[s.fimg]][vdisp_index[s.vdisp]] = s.average[index]
        print s.seqs[0]

    fig = p.figure()
    ax = p3.Axes3D(fig)
    ax.plot_wireframe(x, y, data)
    ax.set_xlabel('Data Weight')
    ax.set_ylabel('Cons. Weight')
    ax.set_zlabel(zlabel)
    print x, y
    
    print 'Data:'
    print data

    p.show()

    

    
    print data
    print fimg, vdisp

f = file(sys.argv[1], 'r');

lines = f.readlines()

stats = []
li = 0
while li < len(lines):
    l = lines[li]
    if l.startswith('-----'):
        result, offset = parse_block(lines[li:])
        stats.append(result)
    li += offset


plot(stats, 0, 'Depth difference')
plot(stats, 1, 'Flow difference')
plot(stats, 2, 'Image difference')
