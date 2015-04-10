


// local header files
#include "../inc/Global.h"
#include "../inc/Makros.h"
#include "../inc/Functions_Shared.h"
#include "../inc/Integrands_pp_ttX.h"

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
  using RunParameters::mScale;
  using RunParameters::mScale2;
  using RunParameters::mt2;
  ////////////////////////////////////////////////////////////////////////////////////////////
  // init ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // create identifier from process id and timestamp
  stringstream run_id;
  run_id << "RUN_" << time(0) << "_" << gettid();

  // create an outpufile
  ofstream log_file;
  if (g_options.logfile) log_file.open(string("results/integration_results_")+=run_id.str());
  // the teestream has problems when no file is openened, so just dump to /dev/null
  else log_file.open("/dev/null");
 
  // create teestream object to tee output to std::cout and the selected file
  teestream couT(std::cout,log_file);
  PRINTS(couT,run_id.str());

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
  // all dimensionful quantities will be normalized to top-mass
  const double MT = 173.34;
  HiggsModel THDM_1;
  // from now on all dimensionful parameters will be normalized to MT
  THDM_1.SetScale(MT);
  // default values
  THDM_1.SetVH(246.0);
  THDM_1.SetMUR(MT);
  THDM_1.SetMUF(MT);
  THDM_1.SetMt(MT);
  THDM_1.SetMb(0.0);
  parse_arguments(argc,argv,g_options,THDM_1);
  // MUR may be changed by user -> set AlphaS after parsing the arguments
  THDM_1.SetAlphaS(pdf_ct10nlo->alphasQ(THDM_1.MUR()*MT));
  THDM_1.Print(couT,MT);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs, etc.
  integrand_par ip;
  ip.higgs_model   = &THDM_1;
  ip.s_hadr        = pow((double)g_options.cme*1000.0,2)/mScale2;
  ip.collect_dist  = g_options.dist;
  ip.pdf           = pdf_ct10nlo; 
  ip.distributions = new HAvec{&MttDistributions,
			       &PT1Distributions,
			       &PT2Distributions,
			       &PT12Distributions,
			       &Y1Distributions,
			       &Y2Distributions,
			       &DYDistributions};
  PRINTS(couT,sqrt(ip.s_hadr)*mScale);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // parameters for VEGAS integration routine
  vegas_par vp;
  vp.calls      = g_options.n_calls;
  vp.verbose    = g_options.verb_level;
  PRINTS(couT,vp.calls);
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // { // set up Higgs boson parameters
  //   using namespace HiggsBosons;
  //   // SetPhi1(g_options.mH,
  //   // 	  g_options.GammaH,
  //   // 	  g_options.At,
  //   // 	  g_options.Bt,
  //   // 	  g_options.Ab,
  //   // 	  g_options.Bb
  //   // 	  );
  //   // // Bernies Higgs [set 2]
  //   // SetPhi1(600.0,
  //   // 	    67.71,
  //   // 	    -1.19,
  //   // 	    1.18,
  //   // 	    0.0,
  //   // 	    0.0
  //   // 	    );
  //   // Bernies Higgs [set 5]
  //   SetPhi1(500.0,
  //   	    59.3,
  //   	    0.0,
  //   	    1.67,
  //   	    0.0,
  //   	    0.0
  //   	    );
  //   SetPhi2(530.0,
  //   	    38.50,//2.245439053359254e+01,
  //   	    1.68,//0.7,
  //   	    0.0,//1.1,
  //   	    0.0,
  //   	    0.0);

  //   PRINTS(couT,Vh*mScale);
  //   PRINTS(couT,M_1*mScale);
  //   PRINTS(couT,G_1*mScale);
  //   PRINTS(couT,At_1*Vh);
  //   PRINTS(couT,Bt_1*Vh);
  //   PRINTS(couT,FH_eff_1/mScale);
  //   PRINTS(couT,FA_eff_1/mScale);
  // }
  ////////////////////////////////////////////////////////////////////////////////
  // { // set up run parameters
  //   using namespace RunParameters;
  //   // set scales
  //   double MU = g_options.ren_scale; // ren./fac. scale in units of mt
  //   SetAlphaS(ip.pdf->alphasQ(MU*mScale));
  //   SetMUR(MU);
  //   SetMUF(MU);
  //   // one or two Higgs bosons?
  //   TwoHDM = 1;
  //   PRINTS(couT,mt*mScale);
  //   PRINTS(couT,AlphaS);
  //   PRINTS(couT,MUR*mScale);
  //   PRINTS(couT,MUF*mScale);
  // }
  ////////////////////////////////////////////////////////////////////////////////
  { // set technical cuts on soft/collinear 2->3 phase space regions
    using namespace Cuts;
    COLL_CUT = 1.0-pow(10.0,-g_options.tech_cut);
    SOFT_CUT = pow(10.0,-g_options.tech_cut);
    PRINTS(couT,SOFT_CUT);
    PRINTS(couT,COLL_CUT);
  }
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  
  // collect integration results and errors
  int const N_res = 6;
  double results[N_res] = {0,0,0,0,0,0};
  double errors [N_res] = {0,0,0,0,0,0};

#ifdef WITH_T_SPIN
  // these are the t and anti-t spin directions in the respective restframes (->only the spatial components relevant)
  // in each event the spin vectors will then be calcualted in the tt-z.m.f and boosted into the labframe (p1,p2-z.m.f.)
  FV S1 = {0.0,1.0,0.0,0.0};
  FV S2 = {0.0,1.0,0.0,0.0};

  couT << endl << "Using spin dependent matrix elements. Spin 4-vectors are:" << endl;
  PRINT_4VEC(S1);
  PRINT_4VEC(S2);
#endif

  Integrator INT(couT);
  
  // // integral over born and virtual
  // Integral Int_2_2_noPDF(2,&IntLimits::Int_lo_2_2[0],&IntLimits::Int_up_2_2[0]);
  // Int_2_2_noPDF.SetIntegrand(&Integrand_2_2);
  // // Bernies sqrt(s_part) limits
  // Int_2_2_noPDF.Lo()[0] = 350.0/mScale;
  // Int_2_2_noPDF.Up()[0] = 800.0/mScale;

  // integral over born and virtual  
  Integral Int_2_2_BV(3,IntLimits::Int_lo_2_2_pdf,IntLimits::Int_up_2_2_pdf);
  Int_2_2_BV.SetIntegrand(&Integrand_2_2_pdf_BV);

  // integral over integrated dipoles 
  Integral Int_2_2_ID(4,IntLimits::Int_lo_2_2_pdf_x,IntLimits::Int_up_2_2_pdf_x);
  Int_2_2_ID.SetIntegrand(&Integrand_2_2_pdf_ID);

  // integral over real - uint. dipoels
  Integral Int_2_3(6,IntLimits::Int_lo_2_3_pdf,IntLimits::Int_up_2_3_pdf);
  Int_2_3.SetIntegrand(&Integrand_2_3_pdf);

  // integral over real - uint. dipoels
  Integral Int_2_3_qg_qq(6,IntLimits::Int_lo_2_3_pdf,IntLimits::Int_up_2_3_pdf);
  Int_2_3_qg_qq.SetIntegrand(&Integrand_2_3_qg_qq_pdf);
  ////////////////////////////////////////
  // this is the upper limit of the M12
  // integration in the 2->3 phase space integral
  // needs to be set accourding to hadronic cme
  Int_2_3.Up()[3]       = sqrt(ip.s_hadr);
  Int_2_3_qg_qq.Up()[3] = sqrt(ip.s_hadr);
  ///////////////////////////////////////



  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 000 00 01  [LO: QCD]
  ////////////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_B_QCD)
    {
      couT << "\n [LO (QCD)]\n";       
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(ip.eval_flags);
      SET_EVAL_B_QCDxQCD(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO QCDxQCD]"));
      
      // activate LO (QCD) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(00 000 1));
      	}

      INT.Integrate(Int_2_2_BV,ip,vp);
      results[0] = vp.result;
      errors[0]  = vp.error;
    }

  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 000 00 10  [LO: PHI+INT]
  ////////////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_B_PHI)
    {
      couT << "\n [LO (PHI+INT)]\n";      
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(ip.eval_flags);
      SET_EVAL_B_PHIxQCD(ip.eval_flags);
      SET_EVAL_B_PHIxPHI(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO PHIxQCD,PHIxPHI]"));
  
      // activate LO (PHI) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(00 001 0));
      	}

      INT.Integrate(Int_2_2_BV,ip,vp);
      results[1] = vp.result;
      errors[1]  = vp.error;
    }

  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 000 01 00  [NLO (V): PHI+INT]
  ////////////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_V)
    {
      couT << "\n [NLO (V)]\n";
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(ip.eval_flags);
      SET_EVAL_V_ALL(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [NLO virt.]"));
 
      // activate NLO (V) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(00 100 0));
      	}

      INT.Integrate(Int_2_2_BV,ip,vp);
      results[2] = vp.result;
      errors[2]  = vp.error;
    }

  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 000 10 00  [NLO (int. Dip.): PHI+INT]
  ////////////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_D)
    {
      couT << "\n [NLO (int. Dip.)]\n";
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(ip.eval_flags);
      SET_EVAL_B_PHIxQCD(ip.eval_flags);
      SET_EVAL_B_PHIxPHI(ip.eval_flags);
      SET_EVAL_D_ALL(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [NLO int. dip.]"));
  
      // activate NLO (ID) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(01 000 0));
      	}

      INT.Integrate(Int_2_2_ID,ip,vp);
      results[3] = vp.result;
      errors[3]  = vp.error;
    }

  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 001 00 00  [NLO (R[GG] - uint. Dip.): PHI+INT]
  ////////////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_R_GG)
    {
      couT << "\n [NLO (R GG)]\n"; 
      // set flags to specify which matrix elements will be evaluated
      ip.eval_flags = 0;//F_EVAL_R_GG;
      // USET_FLAG(F_EVAL_R_PHIxPHI_ISR,ip.eval_flags);
      // USET_FLAG(F_EVAL_R_PHIxPHI_FSR,ip.eval_flags);
      SET_FLAG(F_EVAL_R_ISR_FSR,ip.eval_flags); // distributions ?
      // USET_FLAG(F_EVAL_R_FSR_FSR,ip.eval_flags);
      // USET_FLAG(F_EVAL_R_FSR_ISR,ip.eval_flags);
      // USET_FLAG(F_EVAL_R_FSR_INT,ip.eval_flags);
      ip.SetPS(new PS_2_3(mt2,mt2,0.0,"pp->ttg [NLO real,gg]"));
      
      // activate NLO (R) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(10 000 0));
      	}
      
      INT.Integrate(Int_2_3,ip,vp);
      results[4] = vp.result;
      errors[4]  = vp.error;
    }
  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 110 00 00  [NLO (R[QQ/QG] - uint. Dip.): PHI+INT]
  ////////////////////////////////////////////////////////////////////  
  if (int_flags & (I_FLAGS_R_QQ|I_FLAGS_R_QG))
    {
      couT << "\n [NLO (R QQ & QG)]\n"; 
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(ip.eval_flags);
      SET_EVAL_R_PHIxPHI_QG(ip.eval_flags);
      SET_EVAL_R_PHIxPHI_QQ(ip.eval_flags);
      // SET_EVAL_R_PHIxQCD_QG(ip.eval_flags);
      // SET_EVAL_R_PHIxQCD_QQ(ip.eval_flags);     
      ip.SetPS(new PS_2_3(mt2,mt2,0.0,"pp->ttg [NLO real,qq,qg]"));
      
      // activate NLO (R) histograms 
      for (auto e: *ip.distributions)
      	{
      	  e->SetActive(BOOST_BINARY(10 000 0));
      	}

      INT.Integrate(Int_2_3_qg_qq,ip,vp);
      results[5] = vp.result;
      errors[5]  = vp.error;
    } 

  ////////////////////////////////////////////////////////////
  // PRINT RESULTS
  ////////////////////////////////////////////////////////////
  
  double sum = 0.0;
  double err = 0.0;
  for (int i=0;i<N_res;++i)
    {
      sum += results[i];
      err += pow(errors[i],2);
    }
  err = sqrt(err);
  int prec = (int)log10(sum/err)+3;
  
  couT << endl;
  couT << " ============================================================================================= " << endl;
  couT << " ===== TOTAL X-SECTION ======================================================================= " << endl;
  couT << " ============================================================================================= " << endl;
  couT << "  Born [QCD]             = " << setprecision(prec) << results[0] << " +- " << setprecision(2) << errors[0] << " pb " << endl; 
  couT << "  Born [PHI^2+int.]      = " << setprecision(prec) << results[1] << " +- " << setprecision(2) << errors[1] << " pb " << endl; 
  couT << "  Virtuell               = " << setprecision(prec) << results[2] << " +- " << setprecision(2) << errors[2] << " pb " << endl; 
  couT << "  Dipole int.            = " << setprecision(prec) << results[3] << " +- " << setprecision(2) << errors[3] << " pb " << endl; 
  couT << "  Reell(gg) - Dip. uint. = " << setprecision(prec) << results[4] << " +- " << setprecision(2) << errors[4] << " pb " << endl;
  couT << "  Reell(qp) - Dip. uint. = " << setprecision(prec) << results[5] << " +- " << setprecision(2) << errors[4] << " pb " << endl;
  couT << " --------------------------------------------------------------------------------------------- " << endl;
  couT << "  sum                = " << setprecision(prec) << sum << " +- " << setprecision(2) << err << " pb " << endl; 
  couT << " ============================================================================================= " << endl;
      


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  if (g_options.dist)
    {
      bool PLOT_NLO = int_flags & (I_FLAGS_V|I_FLAGS_D|I_FLAGS_R_GG|I_FLAGS_R_QQ|I_FLAGS_R_QG);
      bool WRITE_TO_FILE = g_options.rootfile;
      TApplication TheApp("MyApp",&argc, argv);
      gROOT->Reset();
      gROOT->SetStyle("Plain");
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleAlign(0);
      
      // create a ROOT file to store the results
      stringstream filename;
      // create identifier from process id and timestamp
      filename << "results/distributions_" << gettid() << "_" << time(0);
      if (int_flags & I_FLAGS_B_QCD ) filename << "_BQCD";
      if (int_flags & I_FLAGS_B_PHI ) filename << "_BPHI";
      if (int_flags & I_FLAGS_V ) filename << "_V";
      if (int_flags & I_FLAGS_D ) filename << "_D";
      if (int_flags & I_FLAGS_R_GG ) filename << "_RGG";
      if (int_flags & I_FLAGS_R_QQ ) filename << "_RQQ";
      if (int_flags & I_FLAGS_R_QG ) filename << "_RQG";        
      filename << ".root";
      TFile* fresults = nullptr;
      if (WRITE_TO_FILE) fresults = new TFile(filename.str().c_str(),"NEW");

      if (ip.distributions != nullptr)
	{
	  for (auto e: *ip.distributions)
	    {
	      // normalize distributions to #iterations, #runs (warm-up does not count)
	      e->Normalize(false);
	    }
	}
      
      double norm[4] = {1.0,//results[0],
			1.0,//results[0]+results[1],
			0.0,
			1.0};//sum};
      // 0 : LO_QCD
      // 1 : LO_PHI
      // 3 : V
      // 4 : D
      // 5 : R

      CanvasPtr c0 = MakeCanvas("Top/Antitop invariant mass distributions",1400,1000);
      DrawDistribution(c0,
		       MttDistributions,
		       &norm[0],
		       string("M_{t#bar{t}} [GeV]"),
		       string("#frac{d#sigma}{dM_{t#bar{t}}} [pb/GeV]"),
		       WRITE_TO_FILE,
		       PLOT_NLO);


      DoubleCanvasPtr cd0 = MakeDoubleCanvas("Top transverse momentum/rapidity distribution");
      DrawTwoDistributions(cd0,
      			   PT1Distributions,
      			   Y1Distributions,
      			   &norm[0],
      			   string("p_{T,t} [GeV]"),string("Y_{t}"),
      			   string("#frac{d#sigma}{dp_{T,t}} [pb/GeV]"),
      			   string("#frac{d#sigma}{dY_{t}} [pb]"),
      			   WRITE_TO_FILE, 
      			   PLOT_NLO);

      CanvasPtr c1 = MakeCanvas("Rapidity difference",1000,1000);
      DrawDistribution(c1,
		       DYDistributions,
		       &norm[0],
		       string("#Delta Y_{t#bar{t}}"),
		       string("#frac{d#sigma}{d #Delta Y_{t#bar{t}}} [pb]"),
		       WRITE_TO_FILE,
		       PLOT_NLO);

      CanvasPtr c2 = MakeCanvas("Top+Antitop transverse momentum distributions",1000,1000);
      DrawDistribution(c2,
		       PT12Distributions,
		       &norm[0],
		       string("p_{T,t+#bar{t}} [GeV]"),
		       string("#frac{d#sigma}{dp_{T,t+#bar{t}}} [pb/GeV]"),
		       WRITE_TO_FILE,
		       PLOT_NLO);      

      // CanvasPtr c0 = MakeCanvas("Partonic cross section",1400,1000);
      // DrawDistribution(c0,
      // 		       SPartDistributions,
      // 		       &norm[0],
      // 		       string("#sqrt{#hat{s}} [GeV]"),
      // 		       string("#hat{#sigma}(#hat{s}) [pb]"),
      // 		       WRITE_TO_FILE,
      // 		       PLOT_NLO);
      
      // DoubleCanvasPtr cd1 = MakeDoubleCanvas("Anti-Top transverse momentum/rapidity distribution");
      // DrawTwoDistributions(cd1,
      // 			   PT2Distributions,
      // 			   Y2Distributions,
      // 			   &norm[0],
      // 			   string("p_{T,#bar{t}} [GeV]"),string("Y_{#bar{t}}"),
      // 			   string("#frac{d#sigma}{dp_{T,#bar{t}}} [pb/GeV]"),
      // 			   string("#frac{d#sigma}{dY_{#bar{t}}} [pb]"),
      // 			   WRITE_TO_FILE, 
      // 			   PLOT_NLO);
 

      if (fresults) fresults->Close();
      TheApp.Run();
    }
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  return 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

