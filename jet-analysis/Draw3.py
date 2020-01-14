#!/usr/bin/env python
from matplotlib.patches import Rectangle
from matplotlib import rc
import numpy as np
import random
import matplotlib.pylab as plt
import matplotlib as mpl
plt.rcParams.update({'figure.max_open_warning': 0})
import argparse
import sys
import yoda

parser = argparse.ArgumentParser(description="Example plotting a YODA histo with Matplotlib")
parser.add_argument('--plot', '-p', type=str, default="13000_no_hadr_MMHT2014lo68cl.yoda")
parser.add_argument('--comp', '-c', type=str, default="900_no_hadr_MMHT2014lo68cl.yoda")
parser.add_argument('--pytqq', '-pytqq', type=str, default="Rivetqq.yoda")
parser.add_argument('--pytgg', '-pytgg', type=str, default="Rivetgg.yoda")
parser.add_argument('--fast', '-f', type=int, default=1)
parser.add_argument('--label1', '-l1', type=str, default="13 000 GeV, Herwig")
parser.add_argument('--label2', '-l2', type=str, default="900 GeV, Herwig")

args = parser.parse_args()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class NumpyHist(object):

    def __init__(self, ao):
        if not isinstance(ao, yoda.AnalysisObject):
            raise Exception("ao argument must be a YODA AnalysisObject; this is a %s" % type(ao))
        ## Get annotations
        self.path = ao.path
        #self.annotations = {aname : ao.annotation(aname) for aname in ao.annotations}
        self.annotations = {}
        for aname in ao.annotations:
            self.annotations[aname] = ao.annotation(aname)
        ## Convert to Scatter and set dimensionality & recarray column names
        s = ao.mkScatter()
        names = ['x', 'y', 'exminus', 'explus', 'eyminus', 'eyplus']
        # TODO: also Scatter1D
        if type(s) is yoda.Scatter2D:
            self.dim = 2
        elif type(s) is yoda.Scatter3D:
            self.dim = 3
            names.insert(2, "z")
            names += ['ezminus', 'ezplus']
        else:
            raise RuntimeError("Whoa! If ao doesn't convert to a 2D or 3D scatter, what is it?!")
        ## Put data points into numpy structure
        dtype = {"names": names, "formats": ["f4" for _ in names]}
        self.data = np.zeros(len(s.points), dtype).view(np.recarray)
        for i, p in enumerate(s.points):
            self.data.x[i] = p.x
            self.data.exminus[i] = p.xErrs[0]
            self.data.explus[i]  = p.xErrs[1]
            self.data.y[i] = p.y
            self.data.eyminus[i] = p.yErrs[0]
            self.data.eyplus[i]  = p.yErrs[1]
            if self.dim > 2:
                self.data.z[i] = p.z
                self.data.ezminus[i] = p.zErrs[0]
                self.data.ezplus[i]  = p.zErrs[1]


    def __len__(self):
        return len(self.x)


    @property
    def xedges_sgl(self):
        return np.append(self.xmin, self.xmax[-1])

    @property
    def xedges_dbl(self):
        edges = np.empty((2*len(self.x),), dtype=self.x.dtype)
        edges[0::2] = self.xmin
        edges[1::2] = self.xmax
        return edges


    @property
    def xmin(self):
        return self.x - self.exminus

    @property
    def xmax(self):
        return self.x + self.explus


    @property
    def ymin(self):
        return self.y - self.eyminus

    @property
    def ymax(self):
        return self.y + self.eyplus


    # TODO: automate more with __get/setattr__?

    @property
    def x(self):
        return self.data.x

    @property
    def exminus(self):
        return self.data.exminus

    @property
    def explus(self):
        return self.data.explus


    @property
    def y(self):
        return self.data.y

    @property
    def eyminus(self):
        return self.data.eyminus

    @property
    def eyplus(self):
        return self.data.eyplus


    # TODO: don't provide these / throw helpful errors if only 2D

    @property
    def z(self):
        return self.data.z

    @property
    def ezminus(self):
        return self.data.ezminus

    @property
    def ezplus(self):
        return self.data.ezplus


    # def __getattr__(self, attr):
    #     "Fall back to the data array for attributes not defined on NumpyHist"
    #     return getattr(self.data, attr)


    def same_binning_as(self, other):
        if self.dim != other.dim:
            return False
        if not (other.x == self.x).all() and \
               (other.exminus == self.exminus).all() and \
               (other.explus == self.explus).all():
            return False
        if self.dim == 2:
            return True
        return (other.y == self.y).all() and \
               (other.eyminus == self.eyminus).all() and \
               (other.eyplus == self.eyplus).all()

def MyPlot(ratio=True, title=None, figsize=(8,6)):
    "Make figure and subplot grid layout"

    #if "plt" not in dir():
    #    mpl, plt = setup_mpl()

    # TODO: Eliminate plt? Requires manual work to set up the backend-specific
    # canvas, but would be better for 'more local' memory management
    fig = plt.figure(figsize=figsize)
    #fig = mpl.figure.Figure(figsize=figsize, tight_layout=True)

    if title:
        fig.suptitle(title, horizontalalignment="center", x=0.5,y=0.95)


    ## Make axes. GridSpec may not be available, in which case fall back ~gracefully
    axmain, axratio = None, None
    if ratio:
        try:
            gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3,1], hspace=0)
            axmain = fig.add_subplot(gs[0])
            axmain.hold(True)
            axratio = fig.add_subplot(gs[1], sharex=axmain)
            axratio.hold(True)
            axratio.axhline(1.0, color="gray") #< Ratio = 1 marker line
        except:
            sys.stderr.write("matplotlib.gridspec not available: falling back to plotting without a ratio\n")
            ratio = False
    if not ratio:
        axmain = fig.add_subplot(1,1,1)
        axmain.hold(True)

    return fig, axmain, axratio

def Myplot_hist_on_axes_1d(axmain, axratio, h, href, h2, href2, h3, href3, default_color="red", default_linestyle="-",href_default_color="blue", href_default_linestyle=".", mylabel = "", href_mylabel = "",logopt=0,myxlabel="",myylabel="",extralabel="",extralabel2="", xmin=-999, xmax=999, ymin = -999,ymax = 999):
    h = mk_numpyhist(h)
    href = mk_numpyhist(href)
    h2 = mk_numpyhist(h2)
    href2 = mk_numpyhist(href2)
    h3 = mk_numpyhist(h3)
    href3 = mk_numpyhist(href3)
    print("h ---------------", h)

    lwidth = 1.0

    print("h.xmin ", h.xmin)
    print("h.xmax ", h.xmax)
    print("h.xmax[-1] ", h.xmax[-1])
    print("np.append(h.y, h.y[-1]) ",np.append(h.y, h.y[-1]))
    artists = axmain.step(np.append(h.xmin, h.xmax[-1]), np.append(h.y, h.y[-1]), label = "13000 GeV MMHT2014lo68cl",where="post", color=default_color, linestyle=default_linestyle, linewidth=lwidth)
    artists = axmain.step(np.append(href.xmin, href.xmax[-1]), np.append(href.y, href.y[-1]), label = "900 GeV MMHT2014lo68cl",where="post", color=href_default_color, linestyle=href_default_linestyle, linewidth=lwidth)
    artists = axmain.step(np.append(h2.xmin, h2.xmax[-1]), np.append(h2.y, h2.y[-1]), label = "13000 GeV CT14lo",where="post", color=default_color, linestyle="-.", linewidth=lwidth)
    artists = axmain.step(np.append(href2.xmin, href2.xmax[-1]), np.append(href2.y, href2.y[-1]), label = "900 GeV CT14lo",where="post", color=href_default_color, linestyle="-.", linewidth=lwidth)
    artists = axmain.step(np.append(h3.xmin, h3.xmax[-1]), np.append(h3.y, h3.y[-1]), label = "13000 GeV NNPDF31\_lo\_as\_0118",where="post", color=default_color, linestyle=":", linewidth=lwidth)
    artists = axmain.step(np.append(href3.xmin, href3.xmax[-1]), np.append(href3.y, href3.y[-1]), label = "900 GeV NNPDF31\_lo\_as\_0118",where="post", color=href_default_color, linestyle=":", linewidth=lwidth)
    #extralabel="hi"
    #extralabel2="hi2"
    axmain.legend(frameon=False,fontsize=18)
    axmain.plot([], [], ' ', label=extralabel)
    axmain.legend(frameon=False,fontsize=18)
    axmain.plot([], [], ' ', label=extralabel2)
    #axmain.legend(frameon=False,prop={fontsize: 18})
    #axmain.legend(extra, "My explanatory text")
    MaxY = []
    MaxY.append(max(h.y))
    MaxY.append(max(href.y))
    if logopt == 0:
        axmain.set_ylim([0.0,max(MaxY)*1.1])
    if (ymax != 999) and (ymin !=-999):
        axmain.set_ylim([ymin,ymax])
    if (xmax != 999) and (xmin !=-999):
        axmain.set_xlim([xmin,xmax])
    else:
        axmain.set_xlim([min(h.xmin),max(h.xmax)])
    #axmain.tick_params(labelbottom=False)
    axmain.set_ylabel(myylabel)
    axmain.set_xlabel(myxlabel)
    
    if logopt == 1:
        axmain.set_yscale("log")
    
    ## Legend entry
    label = h.annotations.get("Parton composition", None)
    if label and artists:
        artists[0].set_label(label)

    return artists

def setup_mpl(engine="MPL", font="TeX Gyre Pagella", fontsize=14, mfont=None, textfigs=True):
    """One-liner matplotlib (mpl) setup.

    By default mpl will be configured with the TeX PGF rendering backend, and a
    Palatino-like font for both text and math contexts, using 'lower-case
    numerals' if supported. Setting the engine to 'TEX' will use standard mpl
    rendering, with calls to LaTeX for axis labels and other text; setting it to
    'MPL' will use the built-in mpl MathText renderer. MPL mode only supports a
    limited set of LaTeX macros and does not render as well as TeX, but it is
    faster and can be used to render to an interactive window.

    The font and mfont optional arguments can be used to choose a different text
    font and math font respectively; if mfont is None, it defaults to the same
    as the text font. The textfigs boolean argument can be set false to disable
    the lower-case/text/old-style numerals and use 'upper-case' numerals
    everywhere. These options do not currently apply to the MPL rendering engine.
    """
    # import matplotlib as mpl
    mpl.rcParams.update({
        "text.usetex" : (engine != "MPL"),
        "font.size"   : int(fontsize),
        "font.family" : "serif", #< TODO: make configurable? auto-detect?
        })

    texpreamble = [r"\usepackage{amsmath,amssymb}", r"\usepackage{mathspec}"]
    mfont = mfont if mfont else font
    fontopts = "[Numbers=OldStyle]" if textfigs else ""
    mfontopts = fontopts.replace("]", ",") + "Scale=MatchUppercase" + "]"
    texpreamble.append( r"\setmainfont{fopts}{{{font}}}".format(fopts=fontopts, font=font) )
    texpreamble.append( r"\setmathsfont(Digits,Latin){fopts}{{{font}}}".format(fopts=mfontopts, font=mfont) )

    if engine == "PGF":
        mpl.use("pgf")
        mpl.rcParams["pgf.preamble"] = texpreamble
    elif engine == "TEX":
        mpl.rcParams["tex.preamble"] = texpreamble

    # TODO: do we need plt?
    from matplotlib import pyplot as plt
    return mpl, plt

def mk_numpyhist(h):
    return h if type(h) is NumpyHist else NumpyHist(h)

aos = yoda.read(args.plot)
bos = yoda.read(args.comp)
cos = yoda.read("13000_no_hadr_CT14lo.yoda")
dos = yoda.read("900_no_hadr_CT14lo.yoda")
eos = yoda.read("13000_no_hadr_03122019_NNPDF31_lo_as_0118.yoda")
fos = yoda.read("900_no_hadr_03122019_NNPDF31_lo_as_0118.yoda")

#ListOfNames = ["GluonAndQuarkFractionPt","GluonMulti","GluonAndQuarkMulti","PartonMulti","GluonFractionPt", "OthersThenGluonAndQuarkMulti","QuarkMulti","PartonFractionPt","QuarkFractionPt","PDGID","OthersThenGluonAndQuarkFractionPt"]

ListOfPDGID = ["PDGID"]
ListOfMulti = ["PartonMulti", "GluonAndQuarkMulti", "GluonMulti", "QuarkMulti", "OthersThenGluonAndQuarkMulti"]
ListOfFractionPt = ["PartonFractionPt", "GluonAndQuarkFractionPt", "GluonFractionPt", "QuarkFractionPt", "OthersThenGluonAndQuarkFractionPt"]
ListOfMultiNames = ["Parton Multiplicity", "Quark and Gluon Multiplicity", "Gluon Multiplicity", "Quark Multiplicity", "Others then Gluon and Quark Multiplicity"]
ListOfFractionPtNames = ["Parton Fraction", "Quark and Gluon Fraction", "Gluon Fraction", "Quark Fraction", "Others then Gluon and Quark Fraction"]


plt.rcParams.update({
    "font.size"   : int(18),
    })

h_ListOfFractionPt13000=[0,0,0,0,0]
h_ListOfMulti13000=[0,0,0,0,0]
h_ListOfPDGID13000 = [0]
    
h_ListOfFractionPt900=[0,0,0,0,0]
h_ListOfMulti900=[0,0,0,0,0]
h_ListOfPDGID900 = [0]

h2_ListOfFractionPt13000=[0,0,0,0,0]
h2_ListOfMulti13000=[0,0,0,0,0]
h2_ListOfPDGID13000 = [0]
    
h2_ListOfFractionPt900=[0,0,0,0,0]
h2_ListOfMulti900=[0,0,0,0,0]
h2_ListOfPDGID900 = [0]

h3_ListOfFractionPt13000=[0,0,0,0,0]
h3_ListOfMulti13000=[0,0,0,0,0]
h3_ListOfPDGID13000 = [0]
    
h3_ListOfFractionPt900=[0,0,0,0,0]
h3_ListOfMulti900=[0,0,0,0,0]
h3_ListOfPDGID900 = [0]

for k in aos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h_ListOfFractionPt13000[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h_ListOfMulti13000[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h_ListOfPDGID13000[0] = k

for k in bos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h_ListOfFractionPt900[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h_ListOfMulti900[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h_ListOfPDGID900[0] = k

for k in cos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h2_ListOfFractionPt13000[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h2_ListOfMulti13000[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h2_ListOfPDGID13000[0] = k

for k in dos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h2_ListOfFractionPt900[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h2_ListOfMulti900[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h2_ListOfPDGID900[0] = k

for k in eos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h3_ListOfFractionPt13000[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h3_ListOfMulti13000[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h3_ListOfPDGID13000[0] = k

for k in fos.values():
    if ("Fast" not in str(k)):
        #print(k)
        for j in range(len(ListOfFractionPt)):
            if '/MC_DIJET_PB/' + ListOfFractionPt[j] == k.path:
                h3_ListOfFractionPt900[j] = k
            if '/MC_DIJET_PB/' + ListOfMulti[j] == k.path:
                h3_ListOfMulti900[j] = k
        if ((ListOfPDGID[0] in k.path)):
            h3_ListOfPDGID900[0] = k

COLORS = ["red", "blue", "magenta", "black", "green"]
LSTYLES = ["-", "-", "-.", ":"]

##NRebin = 30
##
##h_ListOfFractionPt13000[0].rebinBy(NRebin, 10,40)
##h_ListOfFractionPt900[0].rebinBy(NRebin, 10,40)
##h2_ListOfFractionPt13000[0].rebinBy(NRebin, 10,40)
##h2_ListOfFractionPt900[0].rebinBy(NRebin, 10,40)
##h3_ListOfFractionPt13000[0].rebinBy(NRebin, 10,40)
##h3_ListOfFractionPt900[0].rebinBy(NRebin, 10,40)
##
##h_ListOfFractionPt13000[0].rebinBy(4,6,10)
##h_ListOfFractionPt900[0].rebinBy(4,6,10)
##h2_ListOfFractionPt13000[0].rebinBy(4,6,10)
##h2_ListOfFractionPt900[0].rebinBy(4,6,10)
##h3_ListOfFractionPt13000[0].rebinBy(4,6,10)
##h3_ListOfFractionPt900[0].rebinBy(4,6,10)

for i in range(len(ListOfFractionPt)):
    h=h_ListOfFractionPt13000[i]#/h_ListOfFractionPt13000[0]
    href=h_ListOfFractionPt900[i]#/h_ListOfFractionPt900[0]
    h2=h2_ListOfFractionPt13000[i]#/h2_ListOfFractionPt13000[0]
    href2=h2_ListOfFractionPt900[i]#/h2_ListOfFractionPt900[0]
    h3=h3_ListOfFractionPt13000[i]#/h3_ListOfFractionPt13000[0]
    href3=h3_ListOfFractionPt900[i]#/h3_ListOfFractionPt900[0]
    ##if (i != 0):
    ##h.rebinBy(NRebin, 10,40)
    ##href.rebinBy(NRebin, 10,40)
    ##h2.rebinBy(NRebin, 10,40)
    ##href2.rebinBy(NRebin, 10,40)
    ##h3.rebinBy(NRebin, 10,40)
    ##href3.rebinBy(NRebin, 10,40)
    ###
    ##h.rebinBy(4, 6,10)
    ##href.rebinBy(4, 6,10)
    ##h2.rebinBy(4, 6,10)
    ##href2.rebinBy(4, 6,10)
    ##h3.rebinBy(4, 6,10)
    ##href3.rebinBy(4, 6,10)
    #print(type(h))
    h=h/h_ListOfFractionPt13000[0]
    href=href/h_ListOfFractionPt900[0]
    h2=h2/h2_ListOfFractionPt13000[0]
    href2=href2/h2_ListOfFractionPt900[0]
    h3=h3/h3_ListOfFractionPt13000[0]
    href3=href3/h3_ListOfFractionPt900[0]
    fig, axmain, axratio = MyPlot(None, ListOfFractionPtNames[i])
    #print "TYPE: ", type(h)
    firstplot = Myplot_hist_on_axes_1d(axmain, axratio, h, href, h2, href2, h3, href3, COLORS[0], LSTYLES[0], COLORS[1], LSTYLES[1], "13000 GeV","900 GeV",logopt=0,myylabel="Fraction [-]",myxlabel="p$_{\mathrm{T}}$ [GeV]", ymin=0, ymax = 1.5) 
    #def Myplot_hist_on_axes_1d(axmain, axratio, h, href, default_color="red", default_linestyle="-",href_default_color="blue", href_default_linestyle=".", mylabel = "",logopt=0,myxlabel="",myylabel="",extralabel="",extralabel2=""):
    #plot = "second"
    #aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[1], LSTYLES[1], plot,1,0,myxlabel="$\\lambda_0^0$",number=1, extralabel="hi") 
    #plt.show()
    plt.rcParams["savefig.format"] = 'pdf'
    plt.savefig(ListOfFractionPt[i])
    plt.rcParams["savefig.format"] = 'png'
    plt.savefig(ListOfFractionPt[i])
    
for i in range(len(ListOfPDGID)):
    h=h_ListOfPDGID13000[i]
    href=h_ListOfPDGID900[i]
    h2=h2_ListOfPDGID13000[i]
    href2=h2_ListOfPDGID900[i]
    h3=h3_ListOfPDGID13000[i]
    href3=h3_ListOfPDGID900[i]
    fig, axmain, axratio = MyPlot(None, ListOfPDGID[i])
    firstplot = Myplot_hist_on_axes_1d(axmain, axratio, h, href, h2, href2, h3, href3, COLORS[0], LSTYLES[0], COLORS[1], LSTYLES[1], "13000 GeV","900 GeV",logopt=1,myylabel="Events",myxlabel="pdgid of partons",xmin=-20,xmax=84,ymax = 1000) 
    #def Myplot_hist_on_axes_1d(axmain, axratio, h, href, default_color="red", default_linestyle="-",href_default_color="blue", href_default_linestyle=".", mylabel = "",logopt=0,myxlabel="",myylabel="",extralabel="",extralabel2=""):
    #plot = "second"
    #aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[1], LSTYLES[1], plot,1,0,myxlabel="$\\lambda_0^0$",number=1, extralabel="hi")
    #plt.show()
    plt.rcParams["savefig.format"] = 'pdf'
    plt.savefig(ListOfPDGID[i])
    plt.rcParams["savefig.format"] = 'png'
    plt.savefig(ListOfPDGID[i])

for i in range(len(ListOfMulti)):
    h=h_ListOfMulti13000[i]
    href=h_ListOfMulti900[i]
    h2=h2_ListOfMulti13000[i]
    href2=h2_ListOfMulti900[i]
    h3=h3_ListOfMulti13000[i]
    href3=h3_ListOfMulti900[i]

    fig, axmain, axratio = MyPlot(None, ListOfMultiNames[i])
    firstplot = Myplot_hist_on_axes_1d(axmain, axratio, h, href, h2, href2, h3, href3, COLORS[0], LSTYLES[0], COLORS[1], LSTYLES[1], "13000 GeV","900 GeV",logopt=1,myylabel="Events",myxlabel="Number of partons") 
    #def Myplot_hist_on_axes_1d(axmain, axratio, h, href, default_color="red", default_linestyle="-",href_default_color="blue", href_default_linestyle=".", mylabel = "",logopt=0,myxlabel="",myylabel="",extralabel="",extralabel2=""):
    #plot = "second"
    #aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[1], LSTYLES[1], plot,1,0,myxlabel="$\\lambda_0^0$",number=1, extralabel="hi") 
    #plt.show()
    plt.rcParams["savefig.format"] = 'pdf'
    plt.savefig(ListOfMulti[i])
    plt.rcParams["savefig.format"] = 'png'
    plt.savefig(ListOfMulti[i])
