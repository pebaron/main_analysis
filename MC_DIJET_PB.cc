#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Tools/RivetFastJet.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/ATLASConePlugin.hh"
#include "fastjet/CMSIterativeConePlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/D0RunIIConePlugin.hh"
#include "fastjet/TrackJetPlugin.hh"
#include "fastjet/JadePlugin.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include <sstream>
#include <stdio.h> 
#include <math.h>  
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "ROOT.h"
#include "TGraph2D.h"
//#include "fastjet/contrib/ModifiedMassDropTagger.hh"

#include <fastjet/config.h>
#if ! (FASTJET_VERSION_NUMBER >= 30100)
#error "This code requires FastJet >= 3.1.0"
#endif

using namespace fastjet; 

namespace Rivet {

  class MC_DIJET_PB : public Analysis {
  public:

    //DEFAULT_RIVET_ANALYSIS_CTOR(MC_ZINCPB);

    vector<double> JET_AVG_PTMINS;    ///< minimal pt for the avg of the 2 hardest jets pt
    const double JET_MIN_PT_FRACTION; ///< 2nd hardest needs to e at least that fraction of hardest
    const double DELTA_RAP_MAX_DIJET; ///< max rapidity difference between the two jets
    const double LOG_SCALE_MAX;       ///< max value of for the log binning (abs)
    const unsigned int nRADII;        ///< number of radii under consideration
    const unsigned int nQs;           ///< number of scale values (Z ptmin) considered
    const double DELTA_RADII;         ///< radius step size
    const double PARTICLE_RAPMAX;     ///< maximal rapidity allowed for particles
    const double JET_RAPMAX;          ///< maximal rapidity allowed for jets

    MC_DIJET_PB()
      : Analysis("MC_DIJET_PB"),
        JET_MIN_PT_FRACTION(0.8),  // 2nd hardest is at least 0.8 * hardest
        DELTA_RAP_MAX_DIJET(1.0),
        LOG_SCALE_MAX(15.0),
        nRADII(5),
        nQs(5),
        DELTA_RADII(0.2),
        PARTICLE_RAPMAX(2.5),
        JET_RAPMAX(1.5)
    {
      // avg of the two hardest has that minimum
      JET_AVG_PTMINS.push_back( 50.0);
      JET_AVG_PTMINS.push_back(100.0);
      JET_AVG_PTMINS.push_back(200.0);
      JET_AVG_PTMINS.push_back(400.0);
      JET_AVG_PTMINS.push_back(800.0);
    }

    void init() {

      FinalState fs(-PARTICLE_RAPMAX, PARTICLE_RAPMAX, 0.0*GeV);
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      addProjection(jet_input, "JET_INPUT");

      declare(fs, "FS");
      JetDefinition AKT2(ee_genkt_algorithm,0.2,-1.0, WTA_modp_scheme);
      JetDefinition AKT4(ee_genkt_algorithm,0.4,-1.0, WTA_modp_scheme);
      JetDefinition AKT6(ee_genkt_algorithm,0.6,-1.0, WTA_modp_scheme);
      JetDefinition AKT8(ee_genkt_algorithm,0.8,-1.0, WTA_modp_scheme);
      JetDefinition AKT10(ee_genkt_algorithm,1.0,-1.0, WTA_modp_scheme);
      ca_wta_recluster = Recluster(JetDefinition(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme),
                                   false, Recluster::keep_only_hardest);
      declare(FastJets(fs, AKT2), "KT02");
      declare(FastJets(fs, AKT4), "KT04");
      declare(FastJets(fs, AKT6), "KT06");
      declare(FastJets(fs, AKT8), "KT08");
      declare(FastJets(fs, AKT10), "KT10");
      
      
      double Nbins = 250;
      double NbinsMulti = 50;
      double MaxMulti = 250;
      _histFastJets02PtSubLeading = bookHisto1D("FastJets02PtSubLeading", 50, 0, 200);
      _histFastJets02PtLeading = bookHisto1D("FastJets02PtLeading", 50, 0, 200);
      _histFastJets02PtReclust = bookHisto1D("FastJets02PtReclust", 50, 0, 200);
      
      _histFastJets02Pt = bookHisto1D("FastJets02Pt", 50, 0, 200);
      _histFastJets02Mult   = bookHisto1D("FastJets02Mult", 25, 0, 50);
      _histFastJets02E   = bookHisto1D("FastJets02E", 25, 0, 50);
      _histFastJets02Eta    = bookHisto1D("FastJets02Eta", 50, -5, 5);
      _histFastJets02Rapidity    = bookHisto1D("FastJets02Rapidity", 50, -5, 5);
      _histFastJets02Phi    = bookHisto1D("FastJets02Phi", 50, 0, TWOPI);

      _histFastJets02MultLam = bookHisto1D("FastJets02MultLam", NbinsMulti, 0.0, MaxMulti);
      _histFastJets02PtLam = bookHisto1D("FastJets02PtLam", Nbins, 0.0, 1.0);
      _histFastJets02LhaLam = bookHisto1D("FastJets02LhaLam", Nbins, 0.0, 1.0);
      _histFastJets02WidthLam = bookHisto1D("FastJets02WidhtLam", Nbins, 0.0, 1.0);
      _histFastJets02MassLam = bookHisto1D("FastJets02MassLam", Nbins, 0.0, 1.0); 

      _histFastJets04PtSubLeading = bookHisto1D("FastJets04PtSubLeading", 50, 0, 200);
      _histFastJets04PtLeading = bookHisto1D("FastJets04PtLeading", 50, 0, 200);
      _histFastJets04PtReclust = bookHisto1D("FastJets04PtReclust", 50, 0, 200);
      _histFastJets04Pt = bookHisto1D("FastJets04Pt", 50, 0, 200);
      _histFastJets04Mult   = bookHisto1D("FastJets04Mult", 25, 0, 50);
      _histFastJets04E   = bookHisto1D("FastJets04E", 25, 0, 50);
      _histFastJets04Eta    = bookHisto1D("FastJets04Eta", 50, -5, 5);
      _histFastJets04Rapidity    = bookHisto1D("FastJets04Rapidity", 50, -5, 5);
      _histFastJets04Phi    = bookHisto1D("FastJets04Phi", 50, 0, TWOPI);

      _histFastJets04MultLam = bookHisto1D("FastJets04MultLam", NbinsMulti, 0, MaxMulti);
      _histFastJets04PtLam = bookHisto1D("FastJets04PtLam", Nbins, 0.0, 1.0);
      _histFastJets04LhaLam = bookHisto1D("FastJets04LhaLam", Nbins, 0.0, 1.0);
      _histFastJets04WidthLam = bookHisto1D("FastJets04WidhtLam", Nbins, 0.0, 1.0);
      _histFastJets04MassLam = bookHisto1D("FastJets04MassLam", Nbins, 0.0, 1.0);

      _histFastJets06PtSubLeading = bookHisto1D("FastJets06PtSubLeading", 50, 0, 200);
      _histFastJets06PtLeading = bookHisto1D("FastJets06PtLeading", 50, 0, 200);
      _histFastJets06PtReclust = bookHisto1D("FastJets06PtReclust", 50, 0, 200);
      _histFastJets06Pt = bookHisto1D("FastJets06Pt", 50, 0, 200);
      _histFastJets06Mult   = bookHisto1D("FastJets06Mult", 25, 0, 50);
      _histFastJets06E   = bookHisto1D("FastJets06E", 25, 0, 50);
      _histFastJets06Eta    = bookHisto1D("FastJets06Eta", 50, -5, 5);
      _histFastJets06Rapidity    = bookHisto1D("FastJets06Rapidity", 50, -5, 5);
      _histFastJets06Phi    = bookHisto1D("FastJets06Phi", 50, 0, TWOPI);

      _histFastJets06MultLam = bookHisto1D("FastJets06MultLam", NbinsMulti, 0, MaxMulti);
      _histFastJets06PtLam = bookHisto1D("FastJets06PtLam", Nbins, 0.0, 1.0);
      _histFastJets06LhaLam = bookHisto1D("FastJets06LhaLam", Nbins, 0.0, 1.0);
      _histFastJets06WidthLam = bookHisto1D("FastJets06WidhtLam", Nbins, 0.0, 1.0);
      _histFastJets06MassLam = bookHisto1D("FastJets06MassLam", Nbins, 0.0, 1.0);

      _histFastJets08PtSubLeading = bookHisto1D("FastJets08PtSubLeading", 50, 0, 200);
      _histFastJets08PtLeading = bookHisto1D("FastJets08PtLeading", 50, 0, 200);
      _histFastJets08PtReclust = bookHisto1D("FastJets08PtReclust", 50, 0, 200);
      _histFastJets08Pt = bookHisto1D("FastJets08Pt", 50, 0, 200);
      _histFastJets08Mult   = bookHisto1D("FastJets08Mult", 25, 0, 50);
      _histFastJets08E   = bookHisto1D("FastJets08E", 25, 0, 50);
      _histFastJets08Eta    = bookHisto1D("FastJets08Eta", 50, -5, 5);
      _histFastJets08Rapidity    = bookHisto1D("FastJets08Rapidity", 50, -5, 5);
      _histFastJets08Phi    = bookHisto1D("FastJets08Phi", 50, 0, TWOPI);

      _histFastJets08MultLam = bookHisto1D("FastJets08MultLam", NbinsMulti, 0, MaxMulti);
      _histFastJets08PtLam = bookHisto1D("FastJets08PtLam", Nbins, 0.0, 1.0);
      _histFastJets08LhaLam = bookHisto1D("FastJets08LhaLam", Nbins, 0.0, 1.0);
      _histFastJets08WidthLam = bookHisto1D("FastJets08WidhtLam", Nbins, 0.0, 1.0);
      _histFastJets08MassLam = bookHisto1D("FastJets08MassLam", Nbins, 0.0, 1.0);

      _histFastJets10PtSubLeading = bookHisto1D("FastJets10PtSubLeading", 50, 0, 200);
      _histFastJets10PtLeading = bookHisto1D("FastJets10PtLeading", 50, 0, 200);
      _histFastJets10PtReclust = bookHisto1D("FastJets10PtReclust", 50, 0, 200);
      _histFastJets10Pt = bookHisto1D("FastJets10Pt", 50, 0, 200);
      _histFastJets10Mult   = bookHisto1D("FastJets10Mult", 25, 0, 50);
      _histFastJets10E   = bookHisto1D("FastJets10E", 25, 0, 50);
      _histFastJets10Eta    = bookHisto1D("FastJets10Eta", 50, -5, 5);
      _histFastJets10Rapidity    = bookHisto1D("FastJets10Rapidity", 50, -5, 5);
      _histFastJets10Phi    = bookHisto1D("FastJets10Phi", 50, 0, TWOPI);

      _histFastJets10MultLam = bookHisto1D("FastJets10MultLam", NbinsMulti, 0.0, MaxMulti);
      _histFastJets10PtLam = bookHisto1D("FastJets10PtLam", Nbins, 0.0, 1.0);
      _histFastJets10LhaLam = bookHisto1D("FastJets10LhaLam", Nbins, 0.0, 1.0);
      _histFastJets10WidthLam = bookHisto1D("FastJets10WidhtLam", Nbins, 0.0, 1.0);
      _histFastJets10MassLam = bookHisto1D("FastJets10MassLam", Nbins, 0.0, 1.0);

      ArrayOfHist[0][0] = _histFastJets02PtSubLeading;
      ArrayOfHist[0][1] = _histFastJets02PtLeading;
      ArrayOfHist[0][2] = _histFastJets02PtReclust;
      ArrayOfHist[0][3] = _histFastJets02Pt;
      ArrayOfHist[0][4] = _histFastJets02Mult;
      ArrayOfHist[0][5] = _histFastJets02E;
      ArrayOfHist[0][6] = _histFastJets02Eta;
      ArrayOfHist[0][7] = _histFastJets02Rapidity;
      ArrayOfHist[0][8] = _histFastJets02Phi;
      ArrayOfHist[0][9] = _histFastJets02MultLam;
      ArrayOfHist[0][10] = _histFastJets02PtLam;
      ArrayOfHist[0][11] = _histFastJets02LhaLam;
      ArrayOfHist[0][12] = _histFastJets02WidthLam;
      ArrayOfHist[0][13] = _histFastJets02MassLam;
      
      ArrayOfHist[1][0] = _histFastJets04PtSubLeading;
      ArrayOfHist[1][1] = _histFastJets04PtLeading;
      ArrayOfHist[1][2] = _histFastJets04PtReclust;
      ArrayOfHist[1][3] = _histFastJets04Pt;
      ArrayOfHist[1][4] = _histFastJets04Mult;
      ArrayOfHist[1][5] = _histFastJets04E;
      ArrayOfHist[1][6] = _histFastJets04Eta;
      ArrayOfHist[1][7] = _histFastJets04Rapidity;
      ArrayOfHist[1][8] = _histFastJets04Phi;
      ArrayOfHist[1][9] = _histFastJets04MultLam;
      ArrayOfHist[1][10] = _histFastJets04PtLam;
      ArrayOfHist[1][11] = _histFastJets04LhaLam;
      ArrayOfHist[1][12] = _histFastJets04WidthLam;
      ArrayOfHist[1][13] = _histFastJets04MassLam;
      
      ArrayOfHist[2][0] = _histFastJets06PtSubLeading;
      ArrayOfHist[2][1] = _histFastJets06PtLeading;
      ArrayOfHist[2][2] = _histFastJets06PtReclust;
      ArrayOfHist[2][3] = _histFastJets06Pt;
      ArrayOfHist[2][4] = _histFastJets06Mult;
      ArrayOfHist[2][5] = _histFastJets06E;
      ArrayOfHist[2][6] = _histFastJets06Eta;
      ArrayOfHist[2][7] = _histFastJets06Rapidity;
      ArrayOfHist[2][8] = _histFastJets06Phi;
      ArrayOfHist[2][9] = _histFastJets06MultLam;
      ArrayOfHist[2][10] = _histFastJets06PtLam;
      ArrayOfHist[2][11] = _histFastJets06LhaLam;
      ArrayOfHist[2][12] = _histFastJets06WidthLam;
      ArrayOfHist[2][13] = _histFastJets06MassLam;
      
      ArrayOfHist[3][0] = _histFastJets08PtSubLeading;
      ArrayOfHist[3][1] = _histFastJets08PtLeading;
      ArrayOfHist[3][2] = _histFastJets08PtReclust;
      ArrayOfHist[3][3] = _histFastJets08Pt;
      ArrayOfHist[3][4] = _histFastJets08Mult;
      ArrayOfHist[3][5] = _histFastJets08E;
      ArrayOfHist[3][6] = _histFastJets08Eta;
      ArrayOfHist[3][7] = _histFastJets08Rapidity;
      ArrayOfHist[3][8] = _histFastJets08Phi;
      ArrayOfHist[3][9] = _histFastJets08MultLam;
      ArrayOfHist[3][10] = _histFastJets08PtLam;
      ArrayOfHist[3][11] = _histFastJets08LhaLam;
      ArrayOfHist[3][12] = _histFastJets08WidthLam;
      ArrayOfHist[3][13] = _histFastJets08MassLam;
      
      ArrayOfHist[4][0] = _histFastJets10PtSubLeading;
      ArrayOfHist[4][1] = _histFastJets10PtLeading;
      ArrayOfHist[4][2] = _histFastJets10PtReclust;
      ArrayOfHist[4][3] = _histFastJets10Pt;
      ArrayOfHist[4][4] = _histFastJets10Mult;
      ArrayOfHist[4][5] = _histFastJets10E;
      ArrayOfHist[4][6] = _histFastJets10Eta;
      ArrayOfHist[4][7] = _histFastJets10Rapidity;
      ArrayOfHist[4][8] = _histFastJets10Phi;
      ArrayOfHist[4][9] = _histFastJets10MultLam;
      ArrayOfHist[4][10] = _histFastJets10PtLam;
      ArrayOfHist[4][11] = _histFastJets10LhaLam;
      ArrayOfHist[4][12] = _histFastJets10WidthLam;
      ArrayOfHist[4][13] = _histFastJets10MassLam;
      

    }

    void analyze(const Event& event) {
      double PtCut = 20.0;
      const int NumberOfRadiuses = 5;
      const double radius[NumberOfRadiuses] = {0.2,0.4,0.6,0.8,1.0};
      string KTRadius[NumberOfRadiuses] = {"KT02","KT04","KT06","KT08","KT10"};
      Jets jetAr[NumberOfRadiuses];
      const VetoedFinalState &fs = applyProjection<VetoedFinalState>(event,"JET_INPUT"); 
      vector<PseudoJet> particles;
      double Q = 0.0;
      
      particles.reserve(fs.particles().size());
      foreach (const Particle &p, fs.particles()){
        particles.push_back(p.pseudojet());
        Q += p.pT();
      }

      for (size_t i = 0; i < NumberOfRadiuses; ++i) {
      jetAr[i] = apply<FastJets>(event,KTRadius[i]).jetsByPt(PtCut*GeV);
      }

      for (size_t i = 0; i < NumberOfRadiuses; ++i) {
        MSG_DEBUG("Total multiplicity = " << jetAr[i].size());
        ArrayOfHist[i][4]->fill(jetAr[i].size(), event.weight());
      Jet LeadingJet;
      Jet SubLeadingJet;

        size_t k = 0;
        foreach (const Jet& jet_test, jetAr[i]) {
        if (k==0){
                LeadingJet = jet_test;
                }
        if (k==1){
                SubLeadingJet = jet_test;
                if ((SubLeadingJet.pT()/LeadingJet.pT())<0.8){
                //-------------------------------------------------
        size_t j = 0;
        foreach (const Jet& jet, jetAr[i]) {
        
        vector<Jet> JetReclust;
        Jet jet1 = ca_wta_recluster(jet);
        JetReclust.push_back(jet1);
        const double pT = jet.pT();
        const double pTreclust = jet1.pT();
        if (j < 2){
        if (abs(jet.rap()) < 1.5) {
        if (j == 0){
                ArrayOfHist[i][1]->fill(pT/GeV, event.weight());
                }
        if (j == 1){
                ArrayOfHist[i][0]->fill(pT/GeV, event.weight());
                }
        j++;
        ArrayOfHist[i][2]->fill(pTreclust/GeV, event.weight());
        ArrayOfHist[i][3]->fill(pT/GeV, event.weight());
        ArrayOfHist[i][8]->fill(jet.phi(), event.weight());
        ArrayOfHist[i][6]->fill(jet.eta() , event.weight());
        ArrayOfHist[i][5]->fill(jet.E()/GeV , event.weight());
        ArrayOfHist[i][7]->fill(jet.rap() , event.weight());
        double lambdaMult = 0.0;
        double lambdaPt = 0.0;
        double lambdaLha = 0.0;
        double lambdaWidth = 0.0;
        double lambdaMass = 0.0;
            for (const Particle& p : jet.particles()) 
            { 
              Vector3 a = jet.p3();
              Vector3 b = p.p3();
              double alpha = a.get(0)*b.get(0)+a.get(1)*b.get(1)+a.get(2)*b.get(2);
              double anorm = sqrt(a.get(0)*a.get(0)+a.get(1)*a.get(1)+a.get(2)*a.get(2));
              double bnorm = sqrt(b.get(0)*b.get(0)+b.get(1)*b.get(1)+b.get(2)*b.get(2));
              double Acosinus = fabs(acos(alpha/(fabs(anorm)*fabs(bnorm))));
              if (Acosinus > 1.0){
                Acosinus = 1.0;  
                }
              if (fabs(alpha-(fabs(anorm)*fabs(bnorm)))<1e-8){
                Acosinus = 0.0;
                lambdaMult=lambdaMult+1.0;
                lambdaPt=lambdaPt+pow((p.E()/jet.E()),2.0);
                }
                else {
                  lambdaLha=lambdaLha+(p.E()/jet.E())*pow(Acosinus/radius[i],0.5);
                  lambdaMult=lambdaMult+1.0;
                  lambdaPt=lambdaPt+pow((p.E()/jet.E()),2.0);
                  lambdaWidth=lambdaWidth+(p.E()/jet.E())*(Acosinus/radius[i]);
                  lambdaMass=lambdaMass+(p.E()/jet.E())*pow(Acosinus/radius[i],2.0);
                  }
              
            }
                ArrayOfHist[i][9]->fill(lambdaMult,event.weight());
                ArrayOfHist[i][10]->fill(lambdaPt,event.weight());
                ArrayOfHist[i][11]->fill(lambdaLha,event.weight());
                ArrayOfHist[i][12]->fill(lambdaWidth,event.weight());
                ArrayOfHist[i][13]->fill(lambdaMass,event.weight());
          }
          }
          }
                //-------------------------------------------------  
                }
                }
        k++;
        }
        }
     }

    void finalize() {
      normalize({_histFastJets02MultLam,_histFastJets02PtLam,_histFastJets02LhaLam,_histFastJets02WidthLam,_histFastJets02MassLam,_histFastJets02PtReclust,_histFastJets02PtSubLeading,_histFastJets02PtLeading,_histFastJets02Pt, _histFastJets02Mult, _histFastJets02E, _histFastJets02Eta, _histFastJets02Rapidity, _histFastJets02Phi});
      normalize({_histFastJets04MultLam,_histFastJets04PtLam,_histFastJets04LhaLam,_histFastJets04WidthLam,_histFastJets04MassLam,_histFastJets04PtReclust,_histFastJets04PtSubLeading,_histFastJets04PtLeading,_histFastJets04Pt, _histFastJets04Mult, _histFastJets04E, _histFastJets04Eta, _histFastJets04Rapidity, _histFastJets04Phi});
      normalize({_histFastJets06MultLam,_histFastJets06PtLam,_histFastJets06LhaLam,_histFastJets06WidthLam,_histFastJets06MassLam,_histFastJets06PtReclust,_histFastJets06PtSubLeading,_histFastJets06PtLeading,_histFastJets06Pt, _histFastJets06Mult, _histFastJets06E, _histFastJets06Eta, _histFastJets06Rapidity, _histFastJets06Phi});
      normalize({_histFastJets08MultLam,_histFastJets08PtLam,_histFastJets08LhaLam,_histFastJets08WidthLam,_histFastJets08MassLam,_histFastJets08PtReclust,_histFastJets08PtSubLeading,_histFastJets08PtLeading,_histFastJets08Pt, _histFastJets08Mult, _histFastJets08E, _histFastJets08Eta, _histFastJets08Rapidity, _histFastJets08Phi});
      normalize({_histFastJets10MultLam,_histFastJets10PtLam,_histFastJets10LhaLam,_histFastJets10WidthLam,_histFastJets10MassLam,_histFastJets10PtReclust,_histFastJets10PtSubLeading,_histFastJets10PtLeading,_histFastJets10Pt, _histFastJets10Mult, _histFastJets10E, _histFastJets10Eta, _histFastJets10Rapidity, _histFastJets10Phi});
    }


  private:
    Recluster ca_wta_recluster;
    Histo1DPtr _histFastJets02MultLam, _histFastJets02PtLam,_histFastJets02LhaLam,_histFastJets02WidthLam,_histFastJets02MassLam,_histFastJets02PtReclust,_histFastJets02PtSubLeading,_histFastJets02PtLeading,_histFastJets02Pt, _histFastJets02Mult, _histFastJets02E, _histFastJets02Eta, _histFastJets02Rapidity, _histFastJets02Phi;
    Histo1DPtr _histFastJets04MultLam, _histFastJets04PtLam,_histFastJets04LhaLam,_histFastJets04WidthLam,_histFastJets04MassLam,_histFastJets04PtReclust,_histFastJets04PtSubLeading,_histFastJets04PtLeading,_histFastJets04Pt, _histFastJets04Mult, _histFastJets04E, _histFastJets04Eta, _histFastJets04Rapidity, _histFastJets04Phi;
    Histo1DPtr _histFastJets06MultLam, _histFastJets06PtLam,_histFastJets06LhaLam,_histFastJets06WidthLam,_histFastJets06MassLam,_histFastJets06PtReclust,_histFastJets06PtSubLeading,_histFastJets06PtLeading,_histFastJets06Pt, _histFastJets06Mult, _histFastJets06E, _histFastJets06Eta, _histFastJets06Rapidity, _histFastJets06Phi;
    Histo1DPtr _histFastJets08MultLam, _histFastJets08PtLam,_histFastJets08LhaLam,_histFastJets08WidthLam,_histFastJets08MassLam,_histFastJets08PtReclust,_histFastJets08PtSubLeading,_histFastJets08PtLeading,_histFastJets08Pt, _histFastJets08Mult, _histFastJets08E, _histFastJets08Eta, _histFastJets08Rapidity, _histFastJets08Phi;
    Histo1DPtr _histFastJets10MultLam, _histFastJets10PtLam,_histFastJets10LhaLam,_histFastJets10WidthLam,_histFastJets10MassLam,_histFastJets10PtReclust,_histFastJets10PtSubLeading,_histFastJets10PtLeading,_histFastJets10Pt, _histFastJets10Mult, _histFastJets10E, _histFastJets10Eta, _histFastJets10Rapidity, _histFastJets10Phi;
    Histo1DPtr ArrayOfHist[5][14];
  };

  DECLARE_RIVET_PLUGIN(MC_DIJET_PB);

}