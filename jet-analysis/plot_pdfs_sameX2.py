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
parser.add_argument('--q2', '-q2', type=int, default=100)

args = parser.parse_args()


name = args.name
q2 = args.q2

p = lhapdf.mkPDF(name)




## Filling a numpy 2D array with xf(x,Q) samples
xs = [x for x in np.logspace(-3, 0, 100)] # for different xs
xs[0] = 0.0

print xs
gluon_xfs = np.empty([len(xs)])
d_xfs = np.empty([len(xs)])
u_xfs = np.empty([len(xs)])
s_xfs = np.empty([len(xs)])
c_xfs = np.empty([len(xs)])
b_xfs = np.empty([len(xs)])
t_xfs = np.empty([len(xs)])
md_xfs = np.empty([len(xs)])
mu_xfs = np.empty([len(xs)])
ms_xfs = np.empty([len(xs)])
mc_xfs = np.empty([len(xs)])
mb_xfs = np.empty([len(xs)])
mt_xfs = np.empty([len(xs)])

for ix, x in enumerate(xs):
    gluon_xfs[ix] = p.xfxQ(21, x, q2)
    d_xfs[ix] = p.xfxQ(1, x, q2)
    u_xfs[ix] = p.xfxQ(2, x, q2)
    s_xfs[ix] = p.xfxQ(3, x, q2)
    c_xfs[ix] = p.xfxQ(4, x, q2)
    b_xfs[ix] = p.xfxQ(5, x, q2)
    t_xfs[ix] = p.xfxQ(6, x, q2)
    md_xfs[ix] = p.xfxQ(-1, x, q2)
    mu_xfs[ix] = p.xfxQ(-2, x, q2)
    ms_xfs[ix] = p.xfxQ(-3, x, q2)
    mc_xfs[ix] = p.xfxQ(-4, x, q2)
    mb_xfs[ix] = p.xfxQ(-5, x, q2)
    mt_xfs[ix] = p.xfxQ(-6, x, q2)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xs,x*gluon_xfs,'k', label="gluon")
ax.plot(xs,x*d_xfs,'r', label="d quark")
ax.plot(xs,x*u_xfs,'g', label="u quark")
ax.plot(xs,x*s_xfs,'b', label="s quark")
ax.plot(xs,x*c_xfs,'c', label="c quark")
ax.plot(xs,x*b_xfs,'m', label="b quark")
ax.plot(xs,x*t_xfs,'y', label="t quark")
ax.plot(xs,x*md_xfs,'r--', label="anti d quark")
ax.plot(xs,x*mu_xfs,'g--', label="anti u quark")
ax.plot(xs,x*ms_xfs,'b--', label="anti s quark")
ax.plot(xs,x*mc_xfs,'c--', label="anti c quark")
ax.plot(xs,x*mb_xfs,'m--', label="anti b quark")
ax.plot(xs,x*mt_xfs,'y--', label="anti t quark")

ax.set_title('Plot of x $\cdot$ PDF '+name+' at Q$^2$ = '+str(q2)+' GeV')
ax.set_xlabel('x - momentum fraction')
ax.set_ylabel('x$^2\cdot$ PDF')
ax.set_xscale("log")
ax.set_ylim(0,2)

ax.legend()
ax.legend(frameon=False)
#plt.show()
finalname = name+"_Qs_"+str(q2)+"_X2"
plt.rcParams["savefig.format"] = 'pdf'
plt.savefig(finalname)
plt.rcParams["savefig.format"] = 'png'
plt.savefig(finalname)
print " "
print "-------------------------- END OF SCRIPT --------------------------"
print " "