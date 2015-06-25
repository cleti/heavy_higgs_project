


// local header files
#include "../inc/Global.h"
#include "../inc/Makros.h"
#include "../inc/Functions_Shared.h"
#include "../inc/Integrator.h"
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
  if (g_options.logfile)
    {
      log_file.open(g_options.logfile_path+string("integration_pp_ttX_results_")+run_id.str());
      if (!log_file.is_open())
	{
	  cout << endl <<  "Could not open logfile. Path correct? ";
	  PRINT(g_options.logfile_path);
	  exit(1);
	}
    }
  // the teestream has problems when no file is openened, so just dump to /dev/null
  else log_file.open("/dev/null");

  
  // create teestream object to tee output to std::cout and the selected file
  teestream couT(std::cout,log_file);
  PRINTS(couT,run_id.str());
  PRINTS(couT,g_options.useK);

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
  THDM_1.Print(couT,mScale);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs, etc.
  integrand_par ip;
  ip.higgs_model   = &THDM_1;
  ip.s_hadr        = pow((double)g_options.cme*1000.0,2)/mScale2;
  ip.collect_dist  = g_options.dist;
  ip.pdf           = pdf_ct10nlo;
  ip.distributions = new DistVec{
    std::make_shared<Distribution>(&Mtt_Histograms,&OBS_M12),
#ifdef WITH_T_SPIN
    std::make_shared<MeanDistribution>(&Chel_Histograms,&OBS_M12,&OBS_HEL12),
    std::make_shared<MeanDistribution>(&Dopen_Histograms,&OBS_M12,&OBS_D12),
    std::make_shared<MeanDistribution>(&B1_Histograms,&OBS_M12,&OBS_B1),
    std::make_shared<MeanDistribution>(&B2_Histograms,&OBS_M12,&OBS_B2),
    std::make_shared<MeanDistribution>(&OCP1_Histograms,&OBS_M12,&OBS_CP1),
#else
    std::make_shared<Distribution>(&PT1_Histograms,&OBS_PT1),
    std::make_shared<Distribution>(&PT12_Histograms,&OBS_PT12),
    std::make_shared<Distribution>(&Y1_Histograms,&OBS_Y1),
    std::make_shared<Distribution>(&Y2_Histograms,&OBS_Y2),
    std::make_shared<Distribution>(&PT2_Histograms,&OBS_PT2),
    std::make_shared<Distribution>(&DY_Histograms,&OBS_DY12),
#endif
  };
  double CME = sqrt(ip.s_hadr)*mScale;
  PRINTS(couT,CME);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // parameters for VEGAS integration routine
  vegas_par vp;
  // will be adjusted seperately in each integration block
  vp.calls      = g_options.n_calls;
  vp.verbose    = g_options.verb_level;
  PRINTS(couT,vp.calls);
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  { // set technical cuts on soft/collinear 2->3 phase space regions
    using namespace Cuts;
    COLL_CUT = 1.0-pow(10.0,-g_options.tech_cut);
    SOFT_CUT = pow(10.0,-g_options.tech_cut);
    PRINTS(couT,SOFT_CUT);
    PRINTS(couT,1-COLL_CUT);
  }
  ////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // try to open file with QCD NLO data
  if (g_options.qcdfile)
    {
      stringstream filename;
#ifdef WITH_T_SPIN 
      filename << "NLO_QCD_Spin_";
#else
      filename << "NLO_QCD_";
#endif
      filename << std::round(CME/1000) << "TeV_mu" << std::round(THDM_1.MUR()*mScale) << ".dat";
      
      PRINT(filename.str());

      H_IndexMap imap_qcd = {
      	{0,H_LO_QCD},
      	{1,H_NLO_QCD},
	//	{2,H_LO_PHI},
      };
      if (!ReadData(
		    g_options.qcdfile_path+filename.str(),
		    ++(ip.distributions->begin()),
		    ip.distributions->end(),
		    4.0/81.0,
		    imap_qcd
		    ))
	{
	  couT << std::endl << " There have beend errors during data file processing. " << std::endl;
	}
      
      // // if (!file) std::make_shared<std::ifstream>(path.string())
      // std::ifstream file(g_options.qcdfile_path+filename.str());
      // while (!file.is_open())
      // 	{
      // 	  FileBrowser fb(g_options.qcdfile_path);
      // 	  std::cout << std::endl << " Could not open file " << g_options.qcdfile_path << std::endl;
      // 	  std::cout << " Pick new input file:" << std::endl;
      // 	  file.open(fb.browse());
      // 	}

      // H_IndexMap imap1 = {
      // 	{0,H_NLO_QCD},
      // 	{2,H_LO_QCD},
      // };

      // if (!ip.distributions->at(1)->GetHistograms()->ReadTable(file))
      // 	{
      // 	  EXIT(1);
      // 	}
      // if (!ip.distributions->at(2)->GetHistograms()->ReadTable(file))
      // 	{
      // 	  EXIT(1);
      // 	}
       
      // ip.distributions->at(1)->GetHistograms()->ReadFile(H_LO_QCD,
      // 							 g_options.qcdfile_path+filename.str(),
      // 							 0,0);
      // ip.distributions->at(1)->GetHistograms()->ReadFile(H_NLO_QCD,
      // 							 g_options.qcdfile_path+filename.str(),
      // 							 0,1);
      // ip.distributions->at(1)->GetHistograms()->ReadFile(H_LO_PHI,
      // 							 g_options.qcdfile_path+filename.str(),
      // 							 0,2);
      // ip.distributions->at(1)->GetHistograms()->Print(cout);
      // ip.distributions->at(2)->GetHistograms()->Print(cout);
      

      
      //g_options.qcdfile = read_qcd_data(ip.distributions,g_options.qcdfile_path,filename.str());
      // if (!g_options.qcdfile)
      // 	{
      // 	  couT << std::endl << " Error reading QCD data file. " << std::endl;
      // 	}
    }
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  
  // collect integration results and errors
  int const N_res = 6;
  double results[N_res] = {0,0,0,0,0,0};
  double errors [N_res] = {0,0,0,0,0,0};

#ifdef WITH_T_SPIN
  couT << "\n\n =================== Using spin dependent matrix elements =================== \n\n";
  couT << "\n Computing the process pp -> ttbar -> e+ e-  + b-jets, LO top/antitop decay.\n\n";
#else
  couT << "\n Computing the process pp -> ttbar.\n\n";
#endif
  if (g_options.dist) couT << "\n Distributions will be computed.\n\n";
  if (THDM_1.UseK())  couT << "\n Using eff. K-factors to rescale NLO contributions.\n\n";

  Integrator INT(couT);
  
  // // integral over born and virtual
  // Integral Int_2_2_noPDF(2,&IntLimits::Int_lo_2_2[0],&IntLimits::Int_up_2_2[0]);
  // Int_2_2_noPDF.SetIntegrand(&Integrand_2_2);
  // // Bernies sqrt(s_part) limits
  // Int_2_2_noPDF.Lo()[0] = 350.0/mScale;
  // Int_2_2_noPDF.Up()[0] = 800.0/mScale;
  
#ifdef WITH_T_SPIN
  // integral over born and virtual  
  Integral Int_2_2_BV(13,IntLimits::Int_lo_2_6_pdf,IntLimits::Int_up_2_6_pdf);
  Int_2_2_BV.SetIntegrand(&Integrand_2_2_pdf_BV);
#else
  // integral over born and virtual  
  Integral Int_2_2_BV(3,IntLimits::Int_lo_2_2_pdf,IntLimits::Int_up_2_2_pdf);
  Int_2_2_BV.SetIntegrand(&Integrand_2_2_pdf_BV);
#endif
 
  // integral over integrated dipoles 
  Integral Int_2_2_ID(4,IntLimits::Int_lo_2_2_pdf_x,IntLimits::Int_up_2_2_pdf_x);
  Int_2_2_ID.SetIntegrand(&Integrand_2_2_pdf_ID);
  
  // integral over real - uint. dipoles
  Integral Int_2_3(6,IntLimits::Int_lo_2_3_pdf,IntLimits::Int_up_2_3_pdf);
  Int_2_3.SetIntegrand(&Integrand_2_3_pdf);

  // integral over real - uint. dipoles
  Integral Int_2_3_qg_qq(6,IntLimits::Int_lo_2_3_pdf,IntLimits::Int_up_2_3_pdf);
  Int_2_3_qg_qq.SetIntegrand(&Integrand_2_3_qg_qq_pdf);
  ////////////////////////////////////////
  // this is the upper limit of the M12
  // integration in the 2->3 phase space integral
  // needs to be set accourding to hadronic cme
  Int_2_3.Up()[3]       = sqrt(ip.s_hadr); // should be sqrt(..) just for testing!!
  Int_2_3_qg_qq.Up()[3] = sqrt(ip.s_hadr);
  ///////////////////////////////////////



  ////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 000 00 01  [LO: QCD]
  ////////////////////////////////////////////////////////////////////
  if ( (int_flags & I_FLAGS_B_QCD) && !g_options.qcdfile)
    {
      couT << "\n [LO (QCD)]\n";       
      // set flags to specify which matrix elements will be evaluated
      ip.eval_flags = 0;
      SET_EVAL_B_QCDxQCD(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO QCDxQCD]"));
      vp.calls      = g_options.n_calls;
#ifdef WITH_T_SPIN
      ip.ps->set_child(0,new PS_2_3(0,0,0,"t->bl+nu"));
      ip.ps->set_child(1,new PS_2_3(0,0,0,"tbar->bl-nu"));
#endif
      
      // activate LO (QCD) histograms
      for (auto e: *ip.distributions)
	{
	  HistArray* hist = e->GetHistograms();
	  // these spin-observalbes receive no contribution from QCD, only stat. fluctuation
	  // better check for the Avg function here? 
	  if ( hist == &OCP1_Histograms
	       // || hist != &OCP2_Histograms
	       || hist == &B1_Histograms
	       || hist == &B2_Histograms 
	       )
	    {
	      hist->DeactivateAll();
	    }
	  else
	    {
	      hist->SetActive(H_LO_QCD);
	    }
	}
      
      INT.Integrate(Int_2_2_BV,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions) e->GetHistograms()->Normalize();
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
      ip.eval_flags = 0;
      SET_EVAL_B_PHIxQCD(ip.eval_flags);
      SET_EVAL_B_PHIxPHI(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [LO PHIxQCD,PHIxPHI]"));
      vp.calls      = g_options.n_calls/10;
#ifdef WITH_T_SPIN
      ip.ps->set_child(0,new PS_2_3(0,0,0,"t->bl+nu"));
      ip.ps->set_child(1,new PS_2_3(0,0,0,"tbar->bl-nu"));
#endif
      
      // activate LO (PHI) histograms 
      for (auto e: *ip.distributions) e->GetHistograms()->SetActive(H_LO_PHI);

      INT.Integrate(Int_2_2_BV,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions) e->GetHistograms()->Normalize();
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
      ip.eval_flags = 0;
      SET_EVAL_V_ALL(ip.eval_flags);
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [NLO virt.]"));
      vp.calls      = g_options.n_calls/10;
#ifdef WITH_T_SPIN
      ip.ps->set_child(0,new PS_2_3(0,0,0,"t->bl+nu"));
      ip.ps->set_child(1,new PS_2_3(0,0,0,"tbar->bl-nu"));
#endif

      
      // activate NLO (V) histograms 
      for (auto e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_V);

      INT.Integrate(Int_2_2_BV,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions)  e->GetHistograms()->Normalize();
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
      ip.eval_flags = 0;
      SET_EVAL_B_PHIxQCD(ip.eval_flags);
      SET_EVAL_B_PHIxPHI(ip.eval_flags);
      SET_FLAG(F_EVAL_D_GG_ALL,ip.eval_flags);
      SET_FLAG(F_EVAL_D_QG_CONT,ip.eval_flags);
      // need to implement integrated QG dipoles!
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [NLO int. dip.]"));
      vp.calls      = g_options.n_calls*10;
  
      // activate NLO (ID) histograms 
      for (auto e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_ID);

      INT.Integrate(Int_2_2_ID,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions) e->GetHistograms()->Normalize();
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
      ip.eval_flags = F_EVAL_R_GG;
      // USET_FLAG(F_EVAL_R_PHIxPHI_FSR,ip.eval_flags);
      // USET_FLAG(F_EVAL_R_PHIxPHI_ISR,ip.eval_flags);
      ip.SetPS(new PS_2_3(mt2,mt2,0.0,"pp->ttg [NLO real,gg]"));
      vp.calls      = g_options.n_calls;
      
      // activate NLO (R) histograms 
      for (auto e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_R);
      
      INT.Integrate(Int_2_3,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions) e->GetHistograms()->Normalize();
      results[4] = vp.result;
      errors[4]  = vp.error;
    }
  ///////////////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 110 00 00  [NLO (R[QQ/QG] - uint. Dip.): PHI+INT]
  ///////////////////////////////////////////////////////////////////////
  if (int_flags & (I_FLAGS_R_QQ|I_FLAGS_R_QG))
    {
      couT << "\n [NLO (R"; 
      // set flags to specify which matrix elements will be evaluated
      ip.eval_flags = 0;
      if (int_flags & I_FLAGS_R_QG)
	{
	  SET_EVAL_R_PHIxPHI_QG(ip.eval_flags);
	  SET_EVAL_R_PHIxQCD_QG(ip.eval_flags);
	  couT << " QG"; 
	}
      if (int_flags & I_FLAGS_R_QQ)
	{
	  SET_EVAL_R_PHIxPHI_QQ(ip.eval_flags);
	  SET_EVAL_R_PHIxQCD_QQ(ip.eval_flags);
	  couT << " QQ"; 
	}
      couT << ")]\n"; 
      ip.SetPS(new PS_2_3(mt2,mt2,0.0,"pp->ttg [NLO real,qq,qg]"));
      vp.calls      = g_options.n_calls;


      // activate NLO (R) histograms 
      for (auto e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_R);

      INT.Integrate(Int_2_3_qg_qq,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto e: *ip.distributions) e->GetHistograms()->Normalize();
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
      if (std::isnan(results[i]) || std::isinf(results[i]))
	{
	  couT << endl << " results are inf/nan, exiting..." << endl;
	  exit(1);
	}
      err += pow(errors[i],2);
      if (std::isnan(errors[i]) || std::isinf(errors[i]))
	{
	  couT << endl << " results are inf/nan, exiting..." << endl;
	  exit(1);
	}
    }
#ifndef WITH_T_SPIN
  // // this is only for the comparison with P.G.
  // err = sqrt(err)/4;
  // sum /= 4;
#endif
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
      
      // create a ROOT file to store the results
      stringstream rootfile_name;
      // create identifier from process id and timestamp
      rootfile_name << "distributions_" << run_id.str();
      if (int_flags & I_FLAGS_B_QCD ) rootfile_name << "_BQCD";
      if (int_flags & I_FLAGS_B_PHI ) rootfile_name << "_BPHI";
      if (int_flags & I_FLAGS_V ) rootfile_name << "_V";
      if (int_flags & I_FLAGS_D ) rootfile_name << "_D";
      if (int_flags & I_FLAGS_R_GG ) rootfile_name << "_RGG";
      if (int_flags & I_FLAGS_R_QQ ) rootfile_name << "_RQQ";
      if (int_flags & I_FLAGS_R_QG ) rootfile_name << "_RQG";        
      rootfile_name << ".root";
      TFile* file = nullptr;
      if (g_options.rootfile)
	{
	  file = TFile::Open((g_options.rootfile_path+rootfile_name.str()).c_str(),"NEW");
	  // ask user for new path if opening failed, so that nothing gets lost
	  while (!file)
	    {
	      cout << endl << " Could not open ROOT file. Enter new path: " << endl; 
	      cout << " >> "; cin >> g_options.rootfile_path;
	      // attach trailing '/' if not present
	      if (g_options.rootfile_path.back() != '/') g_options.rootfile_path += '/';
	      file = TFile::Open((g_options.rootfile_path+rootfile_name.str()).c_str(),"NEW");
	    }
	}

      
      TApplication TheApp("MyApp",&argc, argv);
      gROOT->Reset();
      gROOT->SetStyle("Plain");
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleAlign(0);
      
      // set these variables accordingly to normalize distributions to total cross section
      double norm[4] = {
	1.0,  // QCD LO
	1.0,  // QCD NLO
	1.0,  // QCD + PHI LO,
	1.0,  // QCD + PHI NLO
      };

      Mtt_Histograms.Print(couT);
      
      CanvasPtr c0 = MakeCanvas("Top/Antitop invariant mass distributions",1000,1000);
      DrawDistribution(c0,
		       Mtt_Histograms,
		       THDM_1,
		       &norm[0],
		       string("M_{t#bar{t}} [GeV]"),
		       string("#frac{d#sigma}{dM_{t#bar{t}}} [pb/GeV]"),
		       g_options.rootfile,
		       PLOT_NLO);

      
#ifdef WITH_T_SPIN
      Dopen_Histograms.Print(couT);
      OCP1_Histograms.Print(couT);
      B1_Histograms.Print(couT);
      Chel_Histograms.Print(couT);
	
 
      CanvasPtr c1 = MakeCanvas("Lepton/Antilepton opening angle",1000,1000);
      DrawDistribution(c1,
      		       Dopen_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("M_{t#bar{t}} [GeV]"),
      		       string("D_{open} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      CanvasPtr c2 = MakeCanvas("CP-odd triple correlation",1000,1000);
      DrawDistribution(c2,
      		       OCP1_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("M_{t#bar{t}} [GeV]"),
      		       string("O_{CP} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      CanvasPtr c31 = MakeCanvas("Top longitudinal polarization",1000,1000);
      DrawDistribution(c31,
      		       B1_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("M_{t#bar{t}} [GeV]"),
      		       string("B_{1} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      CanvasPtr c32 = MakeCanvas("Antiotop longitudinal polarization",1000,1000);
      DrawDistribution(c32,
      		       B2_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("M_{t#bar{t}} [GeV]"),
      		       string("B_{2} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);    
      CanvasPtr c4 = MakeCanvas("Helicity angle distribution",1000,1000);
      DrawDistribution(c4,
      		       Chel_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("M_{t#bar{t}} [GeV]"),
      		       string("C_{hel} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      
#else
      PT1_Histograms.Print(couT);
      PT2_Histograms.Print(couT);
      PT12_Histograms.Print(couT);
      Y1_Histograms.Print(couT);
      Y2_Histograms.Print(couT);     
      DY_Histograms.Print(couT);


      CanvasPtr c11 = MakeCanvas("Top transverse momentum distributions",1000,1000);
      DrawDistribution(c11,
      		       PT1_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("p_{T,t} [GeV]"),
      		       string("#frac{d#sigma}{d p_{T,t}} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      
      CanvasPtr c12 = MakeCanvas("Antitop transverse momentum distributions",1000,1000);
      DrawDistribution(c12,
      		       PT2_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("p_{T,#bar{t}} [GeV]"),
      		       string("#frac{d#sigma}{d p_{T,#\bar{t}}} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);
      
      CanvasPtr c2 = MakeCanvas("Top+Antitop transverse momentum distributions",1000,1000);
      DrawDistribution(c2,
      		       PT12_Histograms,
      		       THDM_1,
      		       &norm[0],
      		       string("|p_{T,t}+p_{T,#bar{t}}| [GeV]"),
      		       string("#frac{d#sigma}{d|p_{T,t}+p_{T,#bar{t}}|} [pb/GeV]"),
      		       g_options.rootfile,
      		       PLOT_NLO);   

      CanvasPtr c31 = MakeCanvas("Top rapidity dsitribution",1000,1000);
      DrawDistribution(c31,
		       Y1_Histograms,
		       THDM_1,
		       &norm[0],
		       string("y_{t}"),
		       string("#frac{d#sigma}{d y_{t}} [pb]"),
		       g_options.rootfile,
		       PLOT_NLO);

      CanvasPtr c32 = MakeCanvas("Antitop rapidity distribution",1000,1000);
      DrawDistribution(c32,
		       Y2_Histograms,
		       THDM_1,
		       &norm[0],
		       string("y_{#bar{t}}"),
		       string("#frac{d#sigma}{d y_{#bar{t}}} [pb]"),
		       g_options.rootfile,
		       PLOT_NLO);
      
      CanvasPtr c4 = MakeCanvas("Top/Antitop rapidity difference distribution",1000,1000);
      DrawDistribution(c4,
		       DY_Histograms,
		       THDM_1,
		       &norm[0],
		       string("#Delta |y|"),
		       string("#frac{d#sigma}{d #Delta |y|} [pb]"),
		       g_options.rootfile,
		       PLOT_NLO);

#endif

      if (file)
	{
	  file->Close();
	  delete file;
	}
      TheApp.Run();
    }
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  return 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

