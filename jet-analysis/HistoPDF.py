#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits import mplot3d
#import matplotlib.pyplot as plt
import lhapdf
import argparse
import yoda
from Draw4_theory_exp import ReturnHisto

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

np.random.seed(19680801)

pt = np.arange(1,201,1,dtype=np.float)
Q2 = 4*pt*pt
x13000 = 2*pt/13000
x900 = 2*pt/900
#print(x900, x13000)
p = lhapdf.mkPDF("MMHT2014nlo68cl")

gluon900 = np.arange(1,201,1,dtype=np.float)
uquark900 = np.arange(1,201,1,dtype=np.float)
dquark900 = np.arange(1,201,1,dtype=np.float)
udquark900 = np.arange(1,201,1,dtype=np.float)
Allquark900 = np.arange(1,201,1,dtype=np.float)

gluon13000 = np.arange(1,201,1,dtype=np.float)
uquark13000 = np.arange(1,201,1,dtype=np.float)
dquark13000 = np.arange(1,201,1,dtype=np.float)
udquark13000 = np.arange(1,201,1,dtype=np.float)
Allquark13000 = np.arange(1,201,1,dtype=np.float)

print("gluon 900")
for i in range(len(pt)):
    #print  pt[i] ,p.xfxQ2(21, x900[i], 	Q2[i])/x900[i]
    gluon900[i] = (p.xfxQ2(21, x900[i], 	Q2[i])/x900[i])
print("u quark 900")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(2, x900[i], 	Q2[i])/x900[i]
    uquark900[i] = (p.xfxQ2(2, x900[i], 	Q2[i])/x900[i])
print("d quark 900")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(1, x900[i], 	Q2[i])/x900[i]
    dquark900[i] = (p.xfxQ2(1, x900[i], 	Q2[i])/x900[i])
print("u+d quark 900")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(1, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(1, x900[i], 	Q2[i])/x900[i] 
    udquark900[i] = (p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] )
print("all quarks 900")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] 
    Allquark900[i] = (p.xfxQ2(1, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(2, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(3, x900[i], 	Q2[i])/x900[i] + p.xfxQ2(-3, x900[i], 	Q2[i])/x900[i]+ p.xfxQ2(-2, x900[i], 	Q2[i])/x900[i]+ p.xfxQ2(-1, x900[i], 	Q2[i])/x900[i] + (p.xfxQ2(21, x900[i], 	Q2[i])/x900[i]))


print("gluon 13000")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(21, x13000[i], 	Q2[i])/x13000[i]
    gluon13000[i] = (p.xfxQ2(21, x13000[i], 	Q2[i])/x13000[i])
print("u quark 13000")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i]
    uquark13000[i] = (p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i])
print("d quark 13000")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(1, x13000[i], 	Q2[i])/x13000[i]
    dquark13000[i] = (p.xfxQ2(1, x13000[i], 	Q2[i])/x13000[i])
print("u+d quark 13000")
for i in range(len(pt)):
    #print pt[i],p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] 
    udquark13000[i] = (p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i] )
print("all quarks 13000")
for i in range(len(pt)):
    Allquark13000[i] = (p.xfxQ2(1, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(2, x13000[i], 	Q2[i])/x13000[i]  + p.xfxQ2(3, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(-3, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(-2, x13000[i], 	Q2[i])/x13000[i] + p.xfxQ2(-1, x13000[i], 	Q2[i])/x13000[i] + (p.xfxQ2(21, x13000[i], 	Q2[i])/x13000[i]))


#Energy = 900
#Energy = 13000 #HERE is the swithcer

##
### Energy 900 GeV
##if Energy == 900:
##    g = np.array([ 0.86 , 0.83 , 0.79, 0.75, 0.72, 0.68, 0.59 , 0.22 ] )
##    u = np.array([0.07,0.09,0.11,0.14,0.16,0.19,0.25,0.53])
##    d = np.array([0.06,0.08,0.09,0.11,0.12,0.13,0.16,0.25])
##
##
### Enregy 13 000 GeV 
##if Energy == 13000:
##    g = np.array([0.89,0.90,0.90,0.89,0.89,0.89,0.88,0.82])
##    u = np.array([0.06,0.05,0.05,0.05,0.06,0.06,0.06,0.10])
##    d = np.array([0.06,0.05,0.05,0.05,0.05,0.05,0.06,0.08])
##
##b = [ 2 , 6 , 10, 14, 18, 22, 32 , 60 ] 
##
fig, ax = plt.subplots(figsize=(15,10))

#ax.plot(pt,gluon900/(gluon900 + udquark900 + udquark900),'--',label='gluon 900 GeV', color='black')
#ax.plot(pt,gluon13000/(gluon13000+udquark13000 + uquark13000),'-',label='gluon 13 000 GeV', color='black')
#ax.plot(pt,(udquark900+uquark900)/(gluon900 + udquark900 + uquark900),'--',label='2u + d quark 900 GeV',color='red')
#ax.plot(pt,(udquark13000+uquark13000)/(gluon13000+udquark13000+ uquark13000),'-',label='2u + d quark 13 000 GeV',color='red')
#ax.plot(b,d,'x',label='d quark',color='red')


#ax.plot(pt,gluon13000/Allquark13000,'-',label='gluon 13 000 GeV', color='black')
#ax.plot(pt,(udquark900+uquark900)/(gluon900 + udquark900 + uquark900),'--',label='quarks 900 GeV',color='red')
#ax.plot(pt,(udquark13000+uquark13000)/(gluon13000+udquark13000+ uquark13000),'-',label='quarks 13 000 GeV',color='red')



h900,h13000 = ReturnHisto()

#ax.hist(h_ListOfFractionPt13000[0],'*',label='gluon measure 13 000 GeV', color='black')
#artists = axmain.step(np.append(h.xmin, h.xmax[-1]), np.append(h.y, h.y[-1]), label = "13000 GeV MMHT2014lo68cl",where="post", color=default_color, linestyle=default_linestyle, linewidth=lwidth)
ax.plot(pt,gluon900/Allquark900,'--',label='gluon 900 GeV', color='black')
artists = ax.step(np.append(h900.xmin, h900.xmax[-1]), np.append(h900.y, h900.y[-1]), label = "Herwig 900 GeV",where="post")
ax.plot(pt,gluon13000/Allquark13000,'-',label='gluon 13 000 GeV', color='black')
artists = ax.step(np.append(h13000.xmin, h13000.xmax[-1]), np.append(h13000.y, h13000.y[-1]),label = "Herwig 13000 GeV",where="post")

ax.set_title('Gluon Fraction PDF and Herwig MHT2014nlo68cl as a function of $p_{T}$')
plt.ylabel('PDF  g / (g + quarks + antiquakrs)')
plt.xlabel('$p_{T}$ [GeV]')


plt.legend(loc='upper right')
fig.tight_layout()

plt.rcParams["savefig.format"] = 'pdf'
plt.savefig("RelativePDF")
plt.rcParams["savefig.format"] = 'png'
plt.savefig("RelativePDF")

plt.show()

