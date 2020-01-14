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
parser.add_argument('--scaleexp', '-s', type=float, default=4)

args = parser.parse_args()


name = args.name
particle = args.particle
p = lhapdf.mkPDF(name)


## Filling a numpy 2D array with xf(x,Q) samples
xs = [x for x in np.logspace(-7, 0, 100)] # for different xs
qs = [q for q in np.logspace(1, args.scaleexp, 100)] # for different scale

xs[0] = 0.0
qs[0] = 0.0
print xs
print qs
gluon_xfs = np.empty([len(xs), len(qs)])
for ix, x in enumerate(xs):
    for iq, q in enumerate(qs):
        gluon_xfs[ix,iq] = p.xfxQ(particle, x, q) # I give q, but it returns q2,  https://lhapdf.hepforge.org/classLHAPDF_1_1PDF.html

#qs = np.log10(qs)
X, Y = np.meshgrid(xs, qs)
print "X" , xs
print "Y" , qs
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, gluon_xfs, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
ax.set_title('Surface Plot of PDF '+name+"_pdgid_"+str(particle))
ax.set_xlabel('x - momentum fraction')
ax.set_ylabel('Qs - Scale')
ax.set_zlabel('PDF')
#ax.set_yscale("log")
ax.invert_xaxis()
#ax.view_init(30,220)
ax.view_init(35, 35)
#fig
#plt.show()
plt.rcParams["savefig.format"] = 'pdf'
plt.savefig(name+"_pdgid_"+str(particle))
plt.rcParams["savefig.format"] = 'png'
plt.savefig(name+"_pdgid_"+str(particle))
print " "
print "-------------------------- END OF SCRIPT --------------------------"
print " "
## Version info, search paths, and metadata
#print(lhapdf.version())
#print(lhapdf.__version__)
#lhapdf.pathsPrepend("/path/to/extra/pdfsets")
#print(lhapdf.paths())
# ...
