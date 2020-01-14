#! /usr/bin/env python2

## Python LHAPDF6 usage example

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits import mplot3d
import numpy as np
#import matplotlib.pyplot as plt
import lhapdf
import argparse

parser = argparse.ArgumentParser(description="Example plotting a YODA histo with Matplotlib")
parser.add_argument('--name', '-n', type=str, default="MMHT2014nlo68cl")
parser.add_argument('--particle', '-p', type=int, default=21)
parser.add_argument('--scaleQ2', '-s', type=float, default=4)
parser.add_argument('--X', '-x', type=float, default=0.01)

args = parser.parse_args()


name = args.name
particle = args.particle
scaleQ2 = args.scaleQ2
x = args.X
p = lhapdf.mkPDF(name)

pt = np.arange(1,201,1,dtype=np.float)
Q2 = 4*pt*pt
x13000 = 2*pt/13000
x900 = 2*pt/900
print(x900, x13000)

print("gluon 900")
for i in range(len(pt)):
    print  pt[i] ,p.xfxQ2(21, x900[i], 	Q2[i])/x900[i]
print("u quark 900")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(2, x900[i], 	Q2[i])/x900[i]
print("d quark 900")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(1, x900[i], 	Q2[i])/x900[i]
print("u+d quark 900")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] 

print("gluon 13000")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(21, x13000[i], 	Q2[i])/x13000[i]
print("u quark 13000")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i]
print("d quark 13000")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(1, x13000[i], 	Q2[i])/x13000[i]
print("u+d quark 13000")
for i in range(len(pt)):
    print pt[i],p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] 

##print("gluons")
##print(p.xfxQ2(21, 0.0044444444, 	16)/0.0044444444)
##print(p.xfxQ2(21, 0.0133333333, 	144)/0.0133333333)
##print(p.xfxQ2(21, 0.0222222222, 	400)/0.0222222222)
##print(p.xfxQ2(21, 0.0311111111, 	784)/0.0311111111)
##print(p.xfxQ2(21, 0.04	       ,  1296)/0.04	       )
##print(p.xfxQ2(21, 0.0488888889, 	1936)/0.0488888889)
##print(p.xfxQ2(21, 0.0711111111, 	4096)/0.0711111111)
##print(p.xfxQ2(21, 0.2222222222, 	40000)/0.2222222222)
##print(p.xfxQ2(21, 0.0003076923, 	16)/0.0003076923)
##print(p.xfxQ2(21, 0.0009230769, 	144)/0.0009230769)
##print(p.xfxQ2(21, 0.0015384615, 	400)/0.0015384615)
##print(p.xfxQ2(21, 0.0021538462, 	784)/0.0021538462)
##print(p.xfxQ2(21, 0.0027692308, 	1296)/0.0027692308)
##print(p.xfxQ2(21, 0.0033846154, 	1936)/0.0033846154)
##print(p.xfxQ2(21, 0.0049230769, 	4096)/0.0049230769)
##print(p.xfxQ2(21, 0.0153846154, 	40000)/0.0153846154)




#print(p.xfxQ2(particle, x, scaleQ2))

