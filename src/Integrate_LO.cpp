


// local header files
#include "../inc/Global.h"
#include "../inc/Makros.h"
#include "../inc/Functions_Shared.h"
#include "../inc/Integrator.h"
#include "../inc/Integrands_pp_ttX.h"
#include "../inc/Cuts.h"

// LHAPDF header files
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"

// ROOT header files
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>
#include <cmath>


using namespace std;

#ifdef WITH_NON_FACT_DIAGRAMS
#include "../ext/LoopTools-2.12/build/clooptools.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
  extern void qlinit_();
#ifdef __cplusplus
}
#endif

  

int& int_flags = g_options.int_flags;
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  ////////////////////////////////////////////////////////////////////////////////////////////
  // init ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // create identifier from process id and timestamp
  stringstream run_id;
  run_id << "RUN_" << time(0) << "_" << gettid();
  program_name = argv[0];

  // create object for model parameters
  HiggsModel THDM_1("generic 2HDM");

  // get parameters from command line arguments
  parse_arguments(argc,argv,g_options,THDM_1);
  
  // create an outpufile
  ofstream log_file;
  log_file.open("/dev/null");
  
  // create teestream object to tee output to std::cout and the selected file
  teestream couT(std::cout,log_file);
  
  SPRINT(couT,run_id.str());
  SPRINT(couT,g_options.useK);

  // create the PDFs
  LHAPDF::PDF* pdf_ct10nlo = LHAPDF::mkPDF("CT10nlo", 0);
  
  // init external libraries
  qlinit_();
#ifdef WITH_NON_FACT_DIAGRAMS
  ltini();
#endif
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // input parameters ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // initialize model parameters
  // some aliases
  double const& mt2     = THDM_1.mt2();
  double const& mScale  = THDM_1.Scale();
  double const& mScale2 = THDM_1.Scale2();
  // extract AlphaS from the PDFs
  if (g_options.as > 0.0)
    {
      THDM_1.SetAlphaS(g_options.as);
    }
  else
    {
      THDM_1.SetAlphaS(pdf_ct10nlo->alphasQ(THDM_1.MUR()*mScale));
    }
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs, etc.
  integrand_par ip;
  ip.higgs_model   = &THDM_1;
  ip.s_hadr        = pow((double)g_options.cme*1000.0,2)/mScale2;
  ip.collect_dist  = true;
  ip.pdf           = pdf_ct10nlo;
  
  ip.distributions = new DistVec{
    std::make_shared<Distribution>(&Mtt_Histograms,&OBS_M12),
  };
  
  double CME = sqrt(ip.s_hadr)*mScale;
  SPRINT(couT,CME);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // parameters for VEGAS integration routine
  vegas_par vp;
  // will be adjusted seperately in each integration block
  vp.calls      = g_options.n_calls;
  vp.verbose    = g_options.verb_level;
  SPRINT(couT,vp.calls);
  ////////////////////////////////////////////////////////////////////////////////

  
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  Integrator INT(couT);

  // integral over 2->2 phase space
  Integral Int_2_2(3,IntLimits::Int_lo_2_2_pdf,IntLimits::Int_up_2_2_pdf);
  Int_2_2.SetIntegrand(&Integrand_2_2_pdf_BV);
  

  ////////////////////////////////////////////////////////////////////
  // [LO: QCD]
  ////////////////////////////////////////////////////////////////////
    {
      couT << "\n [LO (QCD)]\n";       
      // set flags to specify which matrix elements will be evaluated
      ip.eval_flags = 0;
      SET_EVAL_B_QCDxQCD(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO QCDxQCD]"));
      
      // activate LO (QCD) histograms
      for (auto &e: *ip.distributions)
	{
	  e->GetHistograms()->SetActive(H_LO_QCD);
	}
      
      INT.Integrate(Int_2_2,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
    }


  
  ip.eval_flags = 0;
  SET_EVAL_B_PHIxQCD(ip.eval_flags);
  SET_EVAL_B_PHIxPHI(ip.eval_flags);
  ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO PHIxQCD,PHIxPHI]"));

  const int NBosons = 5;
  double masses[NBosons] = {400.0,500.0,600.0,700.0,800.0};
  double widths[NBosons] = { 25.0, 35.0, 45.0, 55.0, 65.0};

  
  ////////////////////////////////////////////////////////////////////
  // [LO: PHI+INT]
  ////////////////////////////////////////////////////////////////////
  for (int i=0; i<1; ++i)
  {
    THDM_1.Print(couT,mScale);
    THDM_1.ClearBosons();
    THDM_1.AddBoson(
		    masses[i],
		    widths[i],
		    1.0,
		    1.0,
		    0.0,
		    0.0);
         
    ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO PHIxQCD,PHIxPHI]"));
      
    // activate LO (PHI) histograms 
    for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_LO_PHI);

    INT.Integrate(Int_2_2,ip,vp);


  }

  // normalize entries of activated histogram to bin width
  for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  // start ROOT environment
  TApplication TheApp("MyApp",&argc, argv);
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleAlign(0);

  // Draw all distributions in DistVec
  for (auto &dist: *ip.distributions)
    {
      HistArray* hist = (dist->GetHistograms());
      hist->Print(couT,g_options.verb_level);
      CanvasPtr c = MakeCanvas(*hist,1000,1000);
      DrawDistribution(c,
		       *hist,
		       THDM_1,
		       g_options.rootfile,
		       true);
    }

  TheApp.Run();
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  return 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

