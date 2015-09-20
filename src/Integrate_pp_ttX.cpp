


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
  THDM_1.Print(couT,mScale);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs, etc.
  integrand_par ip;
  ip.higgs_model   = &THDM_1;
  ip.s_hadr        = pow((double)g_options.cme*1000.0,2)/mScale2;
  ip.collect_dist  = g_options.dist;
  PRINT(g_options.dist);
  ip.pdf           = pdf_ct10nlo;
  
  ip.distributions = new DistVec{
#ifdef WITH_T_SPIN
    std::make_shared<Distribution>(&Mtt_Histograms,&OBS_M12),
    std::make_shared<MeanDistribution>(&Chel_Histograms,&OBS_M12,&OBS_HEL12),
    std::make_shared<MeanDistribution>(&Dopen_Histograms,&OBS_M12,&OBS_D12),
    std::make_shared<MeanDistribution>(&B1_Histograms,&OBS_M12,&OBS_B1),
    // std::make_shared<MeanDistribution>(&B2_Histograms,&OBS_M12,&OBS_B2),
    std::make_shared<MeanDistribution>(&OCP1_Histograms,&OBS_M12,&OBS_CP1),
#else
    // QCD NLO without cuts available
    // std::make_shared<Distribution>(&Mtt_Histograms,&OBS_M12),
    // std::make_shared<Distribution>(&PT1_Histograms,&OBS_PT1),
    // std::make_shared<Distribution>(&PT12_Histograms,&OBS_PT12),
    // std::make_shared<Distribution>(&Y1_Histograms,&OBS_Y1),
    // std::make_shared<Distribution>(&Y2_Histograms,&OBS_Y2),
    
    // QCD NLO with Mtt cut available
    std::make_shared<Distribution>(&Y1_Histograms,&OBS_Y1),
    std::make_shared<Distribution>(&T1_Histograms,&OBS_Theta1),
    std::make_shared<Distribution>(&PT1_Histograms,&OBS_PT1),
#endif
  };

  // at the moment all cuts are evaluated on the parton z.m.f. phase space
  // careful with cuts on observables that are not invariant under boosts in z-direction (rapidity,...)
  ip.cuts = new CutVec{
    std::make_shared<Cut>(&OBS_M12,std::initializer_list<double>{390.0/mScale,540.0/mScale}),
    // std::make_shared<Cut>(&OBS_PT1,std::initializer_list<double>{0.0,100.0/mScale}),
    // std::make_shared<Cut>(&OBS_PT2,std::initializer_list<double>{0.0,100.0/mScale}),
    // rapidity cut does not affect distributions as expected!!!
    // std::make_shared<Cut>(&OBS_Y1,std::initializer_list<double>{0.0,1.0}),
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

  ////////////////////////////////////////////////////////////////////////////////
  { // set technical cuts on soft/collinear 2->3 phase space regions
    using namespace Cuts;
    COLL_CUT = 1.0-pow(10.0,-g_options.tech_cut);
    SOFT_CUT = pow(10.0,-g_options.tech_cut);
    SPRINT(couT,SOFT_CUT);
    SPRINT(couT,1-COLL_CUT);
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
      	{0,H_NLO_QCD},
      	//{1,H_NLO_QCD},
	//{2,H_LO_PHI},
      };
      if (!ReadData(
		    g_options.qcdfile_path+filename.str(),
		    (ip.distributions->begin()),
		    ip.distributions->end(),
		    1.0, // normalization factor, branching ratio 4.0/81.0 not always included
		    imap_qcd
		    ))
	{
	  couT << std::endl << " Error processing data file. " << std::endl;
	}
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
  if (g_options.dist)
    {
      couT << "\n Distributions will be computed:\n";
      for (auto &e: *ip.distributions) couT << " --" << e->GetHistograms()->GetName() << std::endl;
      couT << std::endl;
    }
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
      for (auto &e: *ip.distributions)
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
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
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
      vp.calls      = g_options.n_calls;
#ifdef WITH_T_SPIN
      ip.ps->set_child(0,new PS_2_3(0,0,0,"t->bl+nu"));
      ip.ps->set_child(1,new PS_2_3(0,0,0,"tbar->bl-nu"));
#endif
      
      // activate LO (PHI) histograms 
      for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_LO_PHI);

      INT.Integrate(Int_2_2_BV,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
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
      vp.calls      = g_options.n_calls;
#ifdef WITH_T_SPIN
      ip.ps->set_child(0,new PS_2_3(0,0,0,"t->bl+nu"));
      ip.ps->set_child(1,new PS_2_3(0,0,0,"tbar->bl-nu"));
#endif

      
      // activate NLO (V) histograms 
      for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_V);

      INT.Integrate(Int_2_2_BV,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions)  e->GetHistograms()->Normalize();
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
      SET_FLAG(F_EVAL_D_ALL,ip.eval_flags);
      // need to implement integrated QG dipoles!
      ip.SetPS(new PS_2_2(mt2,mt2,"pp->tt [NLO int. dip.]"));
      vp.calls      = g_options.n_calls*10;
  
      // activate NLO (ID) histograms 
      for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_ID);

      INT.Integrate(Int_2_2_ID,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
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
      ip.eval_flags =  F_EVAL_R_GG;
      // USET_FLAG(F_EVAL_R_PHIxPHI_FSR,ip.eval_flags);
      // USET_FLAG(F_EVAL_R_PHIxPHI_ISR,ip.eval_flags);
      ip.SetPS(new PS_2_3(mt2,mt2,0.0,"pp->ttg [NLO real,gg]"));
      vp.calls      = g_options.n_calls;
      
      // activate NLO (R) histograms 
      for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_R);
      for (auto &e: *ip.distributions) e->GetHistograms()->ToggleActive(H_NLO_PHI_ID);
      
      INT.Integrate(Int_2_3,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
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
      for (auto &e: *ip.distributions) e->GetHistograms()->SetActive(H_NLO_PHI_R);

      INT.Integrate(Int_2_3_qg_qq,ip,vp);

      // normalize entries of activated histogram to bin width
      for (auto &e: *ip.distributions) e->GetHistograms()->Normalize();
      results[5] = vp.result;
      errors[5]  = vp.error;
    } 

  ////////////////////////////////////////////////////////////
  // PRINT RESULTS
  ////////////////////////////////////////////////////////////
  
  double sum  = 0.0;
  double err2 = 0.0;
  for (int i=0;i<N_res;++i)
    {
      sum += results[i];
      if (std::isnan(results[i]) || std::isinf(results[i]))
	{
	  couT << endl << " results are inf/nan, exiting..." << endl;
	  exit(1);
	}
      err2 += pow(errors[i],2);
      if (std::isnan(errors[i]) || std::isinf(errors[i]))
	{
	  couT << endl << " results are inf/nan, exiting..." << endl;
	  exit(1);
	}
    }
  // adjust printed # floating point digits to numerical errors
  int prec = (int)log10(sum/sqrt(err2))+3;
  
  couT << endl;
  couT << " ============================================================================================= " << endl;
  couT << " ===== TOTAL X-SECTION ======================================================================= " << endl;
  couT << " ============================================================================================= " << endl;
  couT << "  Born [QCD]             = " << setprecision(prec) << results[0] << " +- " << setprecision(2) << errors[0] << " pb " << endl; 
  couT << "  Born [PHI^2+int.]      = " << setprecision(prec) << results[1] << " +- " << setprecision(2) << errors[1] << " pb " << endl; 
  couT << "  Virtuell               = " << setprecision(prec) << results[2] << " +- " << setprecision(2) << errors[2] << " pb " << endl; 
  couT << "  Dipole int.            = " << setprecision(prec) << results[3] << " +- " << setprecision(2) << errors[3] << " pb " << endl; 
  couT << "  Reell(gg) - Dip. uint. = " << setprecision(prec) << results[4] << " +- " << setprecision(2) << errors[4] << " pb " << endl;
  couT << "  Reell(qp) - Dip. uint. = " << setprecision(prec) << results[5] << " +- " << setprecision(2) << errors[5] << " pb " << endl;
  couT << " --------------------------------------------------------------------------------------------- " << endl;
  couT << "  sum                    = " << setprecision(prec) << sum << " +- " << setprecision(2) << sqrt(err2) << " pb " << endl; 
  couT << " ============================================================================================= " << endl;
      


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  if (g_options.dist)
    {
      bool PLOT_NLO = int_flags & (I_FLAGS_V|I_FLAGS_D|I_FLAGS_R_GG|I_FLAGS_R_QQ|I_FLAGS_R_QG);
      
      // create unique name for the ROOT file from process id and timestamp
      stringstream rootfile_name;
      rootfile_name << "distributions_" << run_id.str();
      if (int_flags & I_FLAGS_B_QCD ) rootfile_name << "_BQCD";
      if (int_flags & I_FLAGS_B_PHI ) rootfile_name << "_BPHI";
      if (int_flags & I_FLAGS_V )     rootfile_name << "_V";
      if (int_flags & I_FLAGS_D )     rootfile_name << "_D";
      if (int_flags & I_FLAGS_R_GG )  rootfile_name << "_RGG";
      if (int_flags & I_FLAGS_R_QQ )  rootfile_name << "_RQQ";
      if (int_flags & I_FLAGS_R_QG )  rootfile_name << "_RQG";        
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
			   PLOT_NLO);
	}

      if (file)
	{
	  file->Close();
	  delete file;
	}


      TCanvas* c1 = new TCanvas("can_r","can_r",500,500);
      c1->cd(1);
      g_hist_r_wgts.Draw("hist");
      
      TCanvas* c2 = new TCanvas("can_uid","can_uid",500,500);
      c2->cd(1);
      g_hist_uid_wgts.Draw("hist");
      
      TheApp.Run();
    }
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  return 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

