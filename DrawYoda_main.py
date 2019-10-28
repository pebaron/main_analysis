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
parser.add_argument('--plot', '-p', type=str, default="13000.yoda")
parser.add_argument('--comp', '-c', type=str, default="900.yoda")
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
    # fig = mpl.figure.Figure(figsize=figsize, tight_layout=True)

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

def Myplot_hist_on_axes_1d(axmain, axratio, h, href=None, default_color="black", default_linestyle="-", mylabel = "",ylimit=2.5,logopt=0,myxlabel="",myylabel="",number=0,extralabel="",extralabel2=""):
    #if "plt" not in dir():
    #    mpl, plt = setup_mpl()

    h = mk_numpyhist(h)
    href = mk_numpyhist(href)
    # TODO: Split into different plot styles: line/filled/range, step/diag/smooth, ...?

    ## Styles
    default_color = h.annotations.get("Color", default_color)
    marker = h.annotations.get("Marker", h.annotations.get("PolyMarker", None)) # <- make-plots translation
    marker = {"*":"o"}.get(marker, marker) # <- make-plots translation
    mcolor = h.annotations.get("LineColor", default_color)
    errbar = h.annotations.get("ErrorBars", None)
    ecolor = h.annotations.get("ErrorBarsColor", default_color)
    line = h.annotations.get("Line", None)
    lcolor = h.annotations.get("LineColor", default_color)
    lstyle = h.annotations.get("LineStyle", default_linestyle)
    lstyle = {"solid":"-", "dashed":"--", "dotdashed":"-.", "dashdotted":"-.", "dotted":":"}.get(lstyle, lstyle) # <- make-plots translation
    lwidth = 1.0
    msize = 7

    ## If no drawing is enabled, default to a step line
    if not any(h.annotations.get(a) for a in ("Marker", "Line", "ErrorBars")):
        line = "step"

    ## Plotting
    # TODO: Split this into different functions for each kind of data preparation (and smoothing as an extra function?)
    artists = None
    if errbar:
        artists = axmain.errorbar(h.x, h.y, xerr=h.exminus, yerr=h.eyminus, color=ecolor, linestyle="none", linewidth=lwidth, capthick=lwidth) # linestyle="-", marker="o",
    if line == "step":
        artists = axmain.step(np.append(h.xmin, h.xmax[-1]), np.append(h.y, h.y[-1]), label = mylabel,where="post", color=lcolor, linestyle=lstyle, linewidth=lwidth)
        #axmain.legend(frameon=False,prop={fontsize: 18})
        if number==1:
            axmain.plot([], [], ' ', label=extralabel)
            #axmain.legend(frameon=False,prop={fontsize: 18})
        if number==3:
            axmain.plot([], [], ' ', label=extralabel2)
            #axmain.legend(frameon=False,prop={fontsize: 18})
        #axmain.legend(extra, "My explanatory text")
        axmain.set_ylim([0.0001,ylimit])
        #axmain.tick_params(labelbottom=False)
        axmain.set_ylabel("$p(\\lambda)$")
        axmain.set_xlabel(myxlabel)
        if logopt == 1:
            axmain.set_yscale("log")
    elif line == "diag":
        artists = axmain.plot(h.x, h.y, color=lcolor, linestyle=lstyle, linewidth=lwidth)
    elif line == "smooth":
        from scipy.interpolate import spline
        xnew = np.linspace(h.x.min(), h.x.max(), 3*len(h))
        ynew = spline(h.x, h.y, xnew)
        artists = axmain.plot(xnew, ynew, color=lcolor, linestyle=lstyle, linewidth=lwidth)
    if marker:
        artists = axmain.plot(h.x, h.y, marker=marker, markersize=msize, linestyle="none", color=mcolor, markeredgecolor=mcolor)

    ## Legend entry
    label = h.annotations.get("Title", None)
    if label and artists:
        artists[0].set_label(label)

    ## Ratio
    #ratioartists = None
    #if href and h is not href:
    #    # TODO: exclude and specify order via RatioIndex
    #    #assert h.same_binning_as(href)
    #    # TODO: log ratio or #sigma deviation
    #    #print href.y
    #    yratios = h.y/href.y
    #    # TODO: Same styling control as for main plot (with Ratio prefix, default to main plot style)
    #    ## Stepped plot
    #    axratio.set_ylim([0,10])
    #    axratio.set_xlabel(myxlabel)
    #    axratio.set_ylabel("Ratio")
    #    #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    #    ratioartists = axratio.step(href.xedges_sgl, np.append(yratios, yratios[-1]), label = "ratio",where="post", color=lcolor, linestyle=lstyle, linewidth=lwidth)
    #    # TODO: Diag plot
    #    # axratio.plot(href["x"], yratios, color="r", linestyle="--")
    #    # TODO: Smoothed plot

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
#####cos = yoda.read(args.pytqq)
#####dos = yoda.read(args.pytgg)

plt.rcParams.update({
    "font.size"   : int(18),
    })


R=["02","04","06","08","10"]

for loopR in R:
    SaveMe=loopR+"ch.png"
    a1=[]
    for k in aos.values():
        if (("MultLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            a1.append(k)
    for k in bos.values():
        if (("MultLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            a1.append(k)
    b1=[]
    for k in aos.values():
        if (("PtLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            b1.append(k)   
    for k in bos.values():
        if (("PtLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            b1.append(k)   
    c1=[]
    for k in aos.values():
        if (("LhaLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            c1.append(k)
    for k in bos.values():
        if (("LhaLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            c1.append(k)
    d1=[]
    for k in aos.values():
        if (("WidhtLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            d1.append(k)
    for k in bos.values():
        if (("WidhtLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            d1.append(k)
    e1=[]
    for k in aos.values():
        if (("MassLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            e1.append(k)
    for k in bos.values():
        if (("MassLam" in k.path)and(loopR in k.path)and("Ch" not in k.path)):
            e1.append(k)
    a=a1
    b=b1
    c=c1
    d=d1
    e=e1

    Adelta=0
    A0help=[]
    A1help=[]
    j=0
    for a0 in a[0].bins:
        A0help.append(a0.area)
    for a1 in a[1].bins:
        A1help.append(a1.area)
    for i in range(len(A0help)):
        if (A1help[i]+A0help[i])!=0:
            Adelta=Adelta+(A1help[i]-A0help[i])*(A1help[i]-A0help[i])/(A1help[i]+A0help[i])
    Adelta=Adelta/2

    Bdelta=0
    A0help=[]
    A1help=[]
    j=0
    for a0 in b[0].bins:
        A0help.append(a0.area)
    for a1 in b[1].bins:
        A1help.append(a1.area)
    for i in range(len(A0help)):
        if (A1help[i]+A0help[i])!=0:
            Bdelta=Bdelta+(A1help[i]-A0help[i])*(A1help[i]-A0help[i])/(A1help[i]+A0help[i])
    Bdelta=Bdelta/2
    Cdelta=0
    A0help=[]
    A1help=[]
    j=0
    for a0 in c[0].bins:
        A0help.append(a0.area)
    for a1 in c[1].bins:
        A1help.append(a1.area)
    for i in range(len(A0help)):
        if (A1help[i]+A0help[i])!=0:
            Cdelta=Cdelta+(A1help[i]-A0help[i])*(A1help[i]-A0help[i])/(A1help[i]+A0help[i])
    Cdelta=Cdelta/2
    Ddelta=0
    A0help=[]
    A1help=[]
    j=0
    for a0 in d[0].bins:
        A0help.append(a0.area)
    for a1 in d[1].bins:
        A1help.append(a1.area)
    for i in range(len(A0help)):
        if (A1help[i]+A0help[i])!=0:
            Ddelta=Ddelta+(A1help[i]-A0help[i])*(A1help[i]-A0help[i])/(A1help[i]+A0help[i])
    Ddelta=Ddelta/2
    Edelta=0
    A0help=[]
    A1help=[]
    j=0
    for a0 in e[0].bins:
        A0help.append(a0.area)
    for a1 in e[1].bins:
        A1help.append(a1.area)
    for i in range(len(A0help)):
        if (A1help[i]+A0help[i])!=0:
            Edelta=Edelta+(A1help[i]-A0help[i])*(A1help[i]-A0help[i])/(A1help[i]+A0help[i])
    Edelta=Edelta/2
    ############################
    #maxy = [multi,pt,lha,width,mass]
    maxy = [0.1,20,10,10,10]
    if args.fast==1:
        textext="$\\Delta_{Herwig}=\\frac{1}{2}\\int\\frac{(p_q(\\lambda)-p_g(\\lambda))^2}{(p_q(\\lambda)+p_g(\\lambda))}d\\lambda =$ "
        textext2="$\\Delta_{Pythia}=\\frac{1}{2}\\int\\frac{(p_q(\\lambda)-p_g(\\lambda))^2}{(p_q(\\lambda)+p_g(\\lambda))}d\\lambda =$ "
        print a
        print "-----------------"
        
        #COLORS = ["red", "magenta", "blue", "black", "green"]
        COLORS = ["red", "blue", "magenta", "black", "green"]
        #LSTYLES = ["-", "--", "-.", ":"]
        #LSTYLES = ["-", "-.", "-", ":"]
        LSTYLES = ["-", "-", "-.", ":"]
        if "2" in loopR:
            myR="0.2"
        elif "4" in loopR:
            myR="0.4"
        elif "6" in loopR:
            myR="0.6"
        elif "8" in loopR:
            myR="0.8"
        elif "1" in loopR:
            myR="1.0"
        title="Multiplicity, $pp\\rightarrow 2j$, R = "+myR
        href=a[0]
        fig, axmain, axratio = MyPlot(None, title)
        i=0
        
        #some_valid_label = False
        #plot="Herwig"
        for ih, h in enumerate(a):
            #print ih, h.path
            #if i == 2 or i == 3:
            #    plot=args.comp
            #    plot = "R = "+myR+", "+ plot
            #    if "Ch" in h.path:
            #        plot = plot +", charged"
            #else:
            if i==0:
                ##plot=args.plot
                if args.plot == "uu-100-C8.yoda":
                    plot = "quark jets, Herwig before update"
                else:
                    plot = args.label1
            if i==2:
                ##plot=args.comp
                if args.plot == "uu-100-C8.yoda":
                    plot = "gluon jets, Herwig before update"
                else:
                    plot = "gluon jets, Herwig"
            if i==1:
                ##plot=args.comp
                plot = args.label2
            if i==3:
                ##plot=args.comp
                plot = "gluon jets, Pythia"
            

            #if "Ch" n h.path:
            #    plot = plot +", charged"
            aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)], plot,maxy[0],0,myxlabel="$\\lambda_0^0$",number=i, extralabel=textext+str(round(Adelta,2))) 
            i=i+1
        #plt.semilogy(t, 
        if args.plot == "uu-100-C8.yoda":
            plt.savefig("oldLambdaMulti"+SaveMe)
        else:
            plt.savefig("LambdaMulti"+SaveMe)
        title="$p_T^D, pp\\rightarrow 2j$, R = "+myR
        href=b[0]
        fig, axmain, axratio = MyPlot(None, title)
        i=0
        plot="Herwig"
        #some_valid_label = False
        for l in range(len(b)):
            b[l].rebin(6)

        for ih, h in enumerate(b):
            if i==0:
                ##plot=args.plot
                if args.plot == "uu-100-C8.yoda":
                    plot = "quark jets, Herwig before update"
                else:
                    plot = "13 000 GeV, Herwig"
            if i==2:
                ##plot=args.comp
                if args.plot == "uu-100-C8.yoda":
                    plot = "gluon jets, Herwig before update"
                else:
                    plot = "gluon jets, Herwig"
            if i==1:
                ##plot=args.comp
                plot = "900 GeV, Herwig"
            if i==3:
                ##plot=args.comp
                plot = "gluon jets, Pythia"
            
            aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)], plot,maxy[1],0,"$\\lambda_0^2$",number=i, extralabel=textext+str(round(Bdelta,2)))
            i=i+1
        #plt.semilogy(t, 
        if args.plot == "uu-100-C8.yoda":
            plt.savefig("oldLambdaPt"+SaveMe)
        else:
            plt.savefig("LambdaPt"+SaveMe)
        title="LHA, $pp\\rightarrow 2j$, R = "+myR
        href=c[0]
        fig, axmain, axratio = MyPlot(None, title)
        i=0
        plot="Herwig"
        #some_valid_label = False
        for l in range(len(c)):
            c[l].rebin(6)
        for ih, h in enumerate(c):
            if i==0:
                ##plot=args.plot
                if args.plot == "uu-100-C8.yoda":
                    plot = "quark jets, Herwig before update"
                else:
                    plot = "13 000 GeV, Herwig"
            if i==2:
                ##plot=args.comp
                if args.plot == "uu-100-C8.yoda":
                    plot = "gluon jets, Herwig before update"
                else:
                    plot = "gluon jets, Herwig"
            if i==1:
                ##plot=args.comp
                plot = "900 GeV, Herwig"
            if i==3:
                ##plot=args.comp
                plot = "gluon jets, Pythia"
            
            aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)],plot,maxy[2],myxlabel="$\\lambda_{0.5}^1$",number=i, extralabel=textext+str(round(Cdelta,2)))
            i=i+1
        #plt.semilogy(t, 
        if args.plot == "uu-100-C8.yoda":
            plt.savefig("oldLambdaLha"+SaveMe)
        else:
            plt.savefig("LambdaLha"+SaveMe)

        title="Width, $pp\\rightarrow 2j$, R = "+myR
        href=d[0]
        fig, axmain, axratio = MyPlot(None, title)
        i=0
        plot="Herwig"
        for l in range(len(d)):
            d[l].rebin(6)
        #some_valid_label = False
        for ih, h in enumerate(d):
            if i==0:
                ##plot=args.plot
                if args.plot == "uu-100-C8.yoda":
                    plot = "quark jets, Herwig before update"
                else:
                    plot = "13 000 GeV, Herwig"
            if i==2:
                ##plot=args.comp
                if args.plot == "uu-100-C8.yoda":
                    plot = "gluon jets, Herwig before update"
                else:
                    plot = "gluon jets, Herwig"
            if i==1:
                ##plot=args.comp
                plot = "900 GeV, Herwig"
            if i==3:
                ##plot=args.comp
                plot = "gluon jets, Pythia"
            
            aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)],plot,maxy[3],myxlabel="$\\lambda_1^1$",number=i, extralabel=textext+str(round(Ddelta,2)))
            i=i+1
        #plt.semilogy(t, 
        if args.plot == "uu-100-C8.yoda":
            plt.savefig("oldLambdaWidth"+SaveMe)
        else:
            plt.savefig("LambdaWidth"+SaveMe)

        title="Mass, $pp\\rightarrow 2j$, R = "+myR
        href=e[0]
        fig, axmain, axratio = MyPlot(None, title)
        i=0
        plot="Herwig"
        for l in range(len(e)):
            e[l].rebin(6)
        #some_valid_label = False
        for ih, h in enumerate(e):
            if i==0:
                ##plot=args.plot
                if args.plot == "uu-100-C8.yoda":
                    plot = "quark jets, Herwig before update"
                else:
                    plot = "13 000 GeV, Herwig"
            if i==2:
                ##plot=args.comp
                if args.plot == "uu-100-C8.yoda":
                    plot = "gluon jets, Herwig before update"
                else:
                    plot = "gluon jets, Herwig"
            if i==1:
                ##plot=args.comp
                plot = "900 GeV, Herwig"
            if i==3:
                ##plot=args.comp
                plot = "gluon jets, Pythia"
            
            aa = Myplot_hist_on_axes_1d(axmain, axratio, h, href, COLORS[ih % len(COLORS)], LSTYLES[ih % len(LSTYLES)], plot,maxy[4],myxlabel="$\\lambda_2^1$",number=i, extralabel=textext+str(round(Edelta,2)))
            i=i+1
        #plt.semilogy(t, 
        if args.plot == "uu-100-C8.yoda":
            plt.savefig("oldLambdaMass"+SaveMe)
        else:
            plt.savefig("LambdaMass"+SaveMe)

