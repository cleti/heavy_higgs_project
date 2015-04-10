
// cstd lib header files
#include <iostream>
#include <cmath>
#include <bitset>
#include <random>
#include <signal.h>
using namespace std;

// local header files
#include "../inc/Global.h"
#include "../inc/Functions_Shared.h"
#include "../inc/Makros.h"
#include "../inc/Functions_V.h"
#include "../inc/Functions_ID.h"
#include "../inc/Functions_R.h"
#include "../inc/Functions_UID.h"
#include "../inc/Functions_tDecay.h"

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



#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])


double cmp_v_weight(gsl_monte_vegas_state* v_state) {
  if (unlikely(v_state == nullptr)) {
    WARNING("function received nullptr");
    return 0.0;
  }

  // re-calculate VEGAS weight if we have a new event
  // jac = tot. integration volume / (no. of evaluations per grid-bin)
  double jac = v_state->jac; 

  // bin of VEGAS grid where the last function call lies
  int*   bin = v_state->bin; 

  // VEGAS weight = volume of grid-bin
  double   vol = 1.0;
  unsigned dim = v_state->dim;
  for (unsigned j = 0; j < dim; ++j)
    {
      double bin_width;
      unsigned k = bin[j];
      if (k == 0)
	{
	  bin_width = COORD(v_state, 1, j);
	}
      else
	{
	  bin_width = COORD(v_state, k + 1, j) - COORD(v_state, k, j);
	}
      vol *= bin_width;
    }
  return jac*vol;
}


void FillDistributions(std::vector< obs_dist_pair > const & obst_dist_vector,
		       int id,
		       FV const& K1,
		       FV const& K2,
		       FV const& KL1,
		       FV const& KL2,
		       double const & wgt)
{
  // tt distributions
  obs_dist_vector[0].first->Fill(id,obs_M12(K1,K2),wgt);
  obs_dist_vector[1].first->Fill(id,obs_PT(K1)    ,wgt);
  obs_dist_vector[2].first->Fill(id,obs_PT(K2)    ,wgt);
  obs_dist_vector[3].first->Fill(id,obs_Y(K1)     ,wgt);
  obs_dist_vector[4].first->Fill(id,obs_Y(K2)     ,wgt);
  // lepton distributions
  obs_dist_vector[6].first->Fill(id,obs_PT(KL1)     ,wgt);
  obs_dist_vector[7].first->Fill(id,obs_PT(KL2)     ,wgt);
  obs_dist_vector[8].first->Fill(id,obs_PHI(KL1,KL2),wgt);  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int& int_flags = g_options.int_flags;


double Integrand_withPDF_2_2(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  using namespace RunParameters;
  integration_par* ip            = static_cast<integration_par*> (arg);
  PS_2_2* ps                     = dynamic_cast<PS_2_2*> (ip->ps);

  // hadronic c.m.e.
  double s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];
  if (ps->set(sqrt(s_part),x[2]))
    {
      
      ////////////////////////////////////////////////////////////////////////////////////////
      // first the top/anti-top decays because the spin vectors are set in these procedures //
      ////////////////////////////////////////////////////////////////////////////////////////
      double tt_decay = 1.0;
      static PS_2_3* ps_t1 = static_cast<PS_2_3*>(ps->get_child(0));
      static PS_2_3* ps_t2 = static_cast<PS_2_3*>(ps->get_child(1));

      if (unlikely(ps_t1 == nullptr)) {ERROR("'ps_t1': nullptr received");}
      if (unlikely(ps_t2 == nullptr)) {ERROR("'ps_t2': nullptr received");}

      // in the end this sould be the top/antitop spin vectors in the p1+p2 z.m.f. 
      FV& S1 = ps->s1();
      FV& S2 = ps->s2();
		  
      if ( ps_t1->set(mt2,x[4] ,x[5] ,x[6] ,x[7] ,x[8] ,0.0,0.0,0.0) &&
	   ps_t2->set(mt2,x[9] ,x[10],x[11],x[12],x[13],0.0,0.0,0.0)    )
	{
	  ///////////////// decay of top quark ////////////////////////////////////////////////
	  tt_decay *= (ps_t1->get_wgt() * Eval_t_blnu(*ps_t1,S1) / (GammaT));
	  // spin vector in top restframe
	  /////////////////////////////////////////////////////////////////////////////////////

	  ///////////////// decay of anti-top quark ///////////////////////////////////////////
	  tt_decay *= (ps_t2->get_wgt() * Eval_t_blnu(*ps_t2,S2) / (GammaT));
	  S2 *= (-1);// spin vector in antitop restframe
	  /////////////////////////////////////////////////////////////////////////////////////
	  // correcting the phase space factors of 2->2 * ( 2->3 * 2->3 )
	  // to match that of 2->6 phase space
	  tt_decay *= mt2/Pi2; 
	}
      else
	{
	  return 0.0;
	}
      
      FV const& k1 = ps->k1();
      FV const& k2 = ps->k2();
       
      // boost from k1 restframe to tt z.m.f.
      static LT boost_k1_RF;
      boost_k1_RF.set_boost(k1,1);
      boost_k1_RF.apply(S1);
      // boost from k2 restframe to tt z.m.f.
      static LT boost_k2_RF;
      boost_k2_RF.set_boost(k2,1);
      boost_k2_RF.apply(S2);
      
      // FV K1 = {1.0,0.0,0.0,0.0};
      // FV K2 = {1.0,0.0,0.0,0.0};
      // boost_k1_RF.apply(K1);
      // boost_k2_RF.apply(K2);

      // PRINT_4VEC(k1);
      // PRINT_4VEC(k2);
      
      // PRINT_4VEC(K1);
      // PRINT_4VEC(K2);
      // exit(1);
      
      double res_b = 0.0;
      double res_v = 0.0;
      double res_d = 0.0;
      double res_dx= 0.0;
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*s_part);

      // gluon PDFs
      // scale MU must be given in [GeV] !!! 
      double f1 = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2 = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
      
      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      double cf = PREF_GG*f1*f2*CONV_mt2i_pbarn;

      if (g_flags_eval_v)
	{
	  // born matrix elements (GG)
	  res_b = Eval_B(*ps,1)*jac_flux_1*cf;

	  // // born matrix elements (QQ)
	  // if (g_flags_eval_v & F_EVAL_B_QCDxQCD)
	  //   {
	  //     static auto quark_pids     = { 5,  4,  3,  2,  1};	  
	  //     double qq_pdf = 0.0;
	  //     // get the quark/antiquark pdfs
	  //     for (auto i: quark_pids)
	  // 	{
	  // 	  double f1_q  = ip->pdf->xfxQ2(+i, x[0], MUF2*mScale2) / x[0];
	  // 	  double f2_qb = ip->pdf->xfxQ2(-i, x[1], MUF2*mScale2) / x[1];
	  // 	  double f1_qb = ip->pdf->xfxQ2(-i, x[0], MUF2*mScale2) / x[0];
	  // 	  double f2_q  = ip->pdf->xfxQ2(+i, x[1], MUF2*mScale2) / x[1];
	  // 	  qq_pdf += (f1_q*f2_qb+f1_qb*f2_q);
	  // 	}
	  //     // the QCD Born amplitudes are symmetric in scattering angle y
	  //     // therefore the PDF factors can be pulled out
	  //     res_b += Eval_B_QQ(*ps)*jac_flux_1*PREF_QQ*qq_pdf*CONV_mt2i_pbarn;
	  //   }
	  
	  // finite part of the virtual corrections (GG)
	  res_v = Eval_V(*ps)*jac_flux_1*cf;
	}
      
      if (g_flags_eval_id)
	{
	  //////////////////////// I recently changed this, it was at position [1] before /////////////
	  ps->set_x(x[3]); // dont forget to hand down the value of x to the ps object !!!
	  //////////////////////////////////////////////////////////////////////////////////////////////
	  
	  // finite part of the integrated dipoles: delta terms and distribution end-point terms
	  res_d = Eval_ID(*ps)*jac_flux_1*cf;

	  // finite part of the integrated dipoles: x-dependent terms of P and K operators
	  double sx = x[3]*s_part;
	  if ( sx > 4.0*mt2 )
	    {
	      // boosted phase space (boost is done in Eval_ID_X)
	      PS_2_2 ps_x = *ps;
	      // modified version of flux and phase space density
	      res_dx = Eval_ID_X(*ps,ps_x);
	      // the ps_x object is modified in Eval_ID_X -> use the pahsespace weight only after applying this function !!!
	      double jac_flux_x = TwoPi*ps_x.get_wgt()/(2.0*sx);
	      res_dx *= jac_flux_x*cf;
	    }
	}



      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (g_options.dist)
	{
	  vector<obs_dist_pair>* dist = static_cast<vector<obs_dist_pair>*>(ip->distributions);
	  if (unlikely(dist==nullptr))
	    {
	      ERROR("You want distributions? Then pass something more than a nullptr.");
	    }

	  FV const& kl1 = ps_t1->k1();
	  FV const& kl2 = ps_t2->k1();
	  
	  double vwgt = cmp_v_weight(ip->v_state);
	  if (g_flags_eval_v & F_EVAL_B_QCDxQCD)
	    { // only QCD LO: fill in the 0-th hostogram
	      FillDistributions(*dist,0,k1,k2,kl1,kl2,res_b*vwgt*tt_decay);
	      // if (dist->at(0).first->IsActive(0))
	      // 	{
	      // 	  rap_plot->Fill(obs_Y(K1),obs_Y(K2),res_b*vwgt*tt_decay);
	      // 	  pdf_plot->Fill(x[0],x[1]);
	      // 	}
	    }
	  else
	    { // PHI + INT LO: fill in the first histogram
	      FillDistributions(*dist,1,k1,k2,kl1,kl2,res_b*vwgt*tt_decay);
	    }
	  FillDistributions(*dist,3,k1,k2,kl1,kl2,res_v*vwgt*tt_decay);
	  FillDistributions(*dist,4,k1,k2,kl1,kl2,res_d*vwgt*tt_decay);
	  FillDistributions(*dist,4,k1,k2,kl1,kl2,res_dx*vwgt*tt_decay);
	}
      ////////////////////////////////////////////////////////////////////////////
      
      return  (res_b+res_v+res_d+res_dx)*tt_decay;
    }
  else
    {
      return 0.0;
    }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double COLL_CUT = 1.0-1e-7;
double SOFT_CUT = 1e-7;

double Integrand_withPDF_2_3(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  using namespace RunParameters;
  static const int gluon_pid = 21;
  // static auto quark_pids = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5};
  
  integration_par* ip            = static_cast<integration_par*> (arg);
  PS_2_3* ps                     = dynamic_cast<PS_2_3*> (ip->ps);
  gsl_monte_vegas_state* v_state = ip->v_state;
  vector< obs_dist_pair>* dist   = static_cast<vector< obs_dist_pair>*>(ip->distributions);


  // hadronic c.m.e.
  double s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];

  // make a technical cut on the angle between p1/p2 and p3 
  if ( ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5]))// && fabs(x[2])<COLL_CUT )
    {
      FV const& p3 = ps->p3();
      // if (p3[0] < SOFT_CUT) return 0.0;


      // Peters implementation of the cut
      FV const& p1 = ps->p1();
      FV const& p2 = ps->p2();
      FV const& k1 = ps->k1();
      FV const& k2 = ps->k2();
      double S31 = sp(k1,p3)/s_part;
      double S32 = sp(k2,p3)/s_part;
      double s31 = sp(p1,p3)/s_part;
      double s32 = sp(p2,p3)/s_part;
      if ( (s31) < SOFT_CUT || (s32) < SOFT_CUT ||  (S31) < SOFT_CUT || (S32) < SOFT_CUT ) return 0;
      
      // // top anti-top 4-vectors given in the tt z.m.f.
      // FV const& k1 = ps->k1();
      // FV const& k2 = ps->k2();
      FV& S1 = ps->s1();
      FV& S2 = ps->s2();
      S1 = ps->s1_r();
      S2 = ps->s2_r();
      
      // boost from k1 restframe to tt z.m.f.
      static LT boost_k1_RF;
      boost_k1_RF.set_boost(k1,1);
      boost_k1_RF.apply(S1);
      // boost from k2 restframe to tt z.m.f.
      static LT boost_k2_RF;
      boost_k2_RF.set_boost(k2,1);
      boost_k2_RF.apply(S2);


      
      double vwgt = (g_options.dist)?(cmp_v_weight(v_state)):(0.0);
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);

      // the PDF wants dimensionful quantities in units of GeV (here, MUF is in units of mt!) 
      double f1 = ip->pdf->xfxQ2(gluon_pid, x[0], MUF2*mScale2) / x[0];
      double f2 = ip->pdf->xfxQ2(gluon_pid, x[1], MUF2*mScale2) / x[1];
      double cf = PREF_GG*f1*f2*CONV_mt2i_pbarn;

      double res_r_gg = Eval_R_gg(*ps)*cf*jac_flux;
      double res_d_gg = Eval_UID (*ps,-vwgt*cf*jac_flux,dist)*cf*jac_flux;
      
      if (dist && (g_options.dist) ) 
      	{
      	  for (auto e: *dist)
      	    {
      	      e.first->Fill(5,e.second(),res_r_gg*vwgt);
      	    }
      	}
      return (res_r_gg-res_d_gg);
    }
  else
    {
      return 0.0;
    }
}




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  using namespace Constants;
  using namespace RunParameters;
  using namespace Bosons;
  using namespace Prefactors;


  ////////////////////////////////////////////////////////////////////////////////////////////
  // init ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  // create identifier from process id and timestamp
  stringstream run_id;
  run_id << "RUN_" << time(0) << "_" << gettid();

  // create an outpufile
  ofstream log_file;
  log_file.open(string("results/integration_results_")+=run_id.str());

  // create teestream object to tee output to std::couT and file
  teestream couT(std::cout,log_file);
  PRINTS(couT,run_id.str());
  
  // init external libraries
  qlinit_();
#ifdef WITH_NON_FACT_DIAGRAMS
  ltini();
#endif

  ////////////////////////////////////////////////////////////////////////////////////////////
  // input parameters ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // parse arguments and initialize input parameters
  parse_arguments(argc,argv,g_options);

  // technical cuts on soft/collinear 2->3 phase space regions
  COLL_CUT = 1.0-pow(10.0,-g_options.tech_cut);
  SOFT_CUT = pow(10.0,-g_options.tech_cut);


  // scalar boson mass, width and couplings
  SetPhi1(g_options.mH,
	  g_options.GammaH,
	  g_options.At,
	  g_options.Bt,
	  g_options.Ab,
	  g_options.Bb
	  );
  SetPhi2(530.0,
	  38.50,//2.245439053359254e+01,
	  1.68,//0.7,
	  0.0,//1.1,
	  0.0,
	  0.0);
  RunParameters::TwoHDM = 0;
  
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs
  integration_par ip; // 
  ip.s_hadr = pow((double)g_options.cme*1000.0,2)/mScale2; // center-of-mass energy
  double const& S_had = ip.s_hadr;
  ip.pdf    = LHAPDF::mkPDF("CT10nlo", 0); 
  ip.distributions = &obs_dist_vector;

  
  double MU = g_options.ren_scale; // ren./fac. scale in units of mt
  SetAlphaS(ip.pdf->alphasQ(MU*mScale));
  SetMUR(MU);
  SetMUF(MU);


  // parameters for VEGAS integration routine
  vegas_par vp = {
    1, // verbose
    g_options.n_calls, // #calls from argument list
    5, // iterations
    10, // max_runs
    0, // num_runs
    0.2, // chisq_limit
    1, // do_warmup
    0 // grid_fixed
  };

  // print parameters
  PRINTS(couT,mt*mScale);
  PRINTS(couT,Vh*mScale);
  PRINTS(couT,AlphaS);
  PRINTS(couT,M_1*mScale);
  PRINTS(couT,G_1*mScale);
  PRINTS(couT,At_1*Vh);
  PRINTS(couT,Bt_1*Vh);
  PRINTS(couT,FH_1/mScale);
  PRINTS(couT,FA_1/mScale);
  PRINTS(couT,MUR*mScale);
  PRINTS(couT,MUF*mScale);
  PRINTS(couT,sqrt(S_had)*mScale);
  PRINTS(couT,SOFT_CUT);
  PRINTS(couT,COLL_CUT);
  PRINTS(couT,vp.calls);
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  // collect integration results and errors
  double results[5] = {0,0,0,0,0};
  double errors[5]  = {0,0,0,0,0};

  ////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 00001  [LO: QCD]
  ////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_B_QCD)
    {
      PS_2_2 ps_gg_tt("gg->tt-bar");
      PS_2_3 ps_t1("decay_ps_t");
      PS_2_3 ps_t2("decay_ps_tb");
      ps_gg_tt.set_child(0,&ps_t1);
      ps_gg_tt.set_child(1,&ps_t2);
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(g_flags_eval_v);
      USET_EVAL_ALL(g_flags_eval_id);
      SET_EVAL_B_QCDxQCD(g_flags_eval_v);

      ip.ps     = &ps_gg_tt;
      vp.calls  = g_options.n_calls;
      // integrate LO matrix elements with PDFs
      gsl_monte_vegas_state* v_state = nullptr;
      
      // only LO histograms get filled
      for (auto e: obs_dist_vector)
	{
	  e.first->SetActive(BOOST_BINARY(00 000 1));
	}

      couT << endl;
      couT << endl << " EVAL_V_FLAGS = " << bitset<17>(g_flags_eval_v).to_string() << endl;
      couT << endl << " EVAL_D_FLAGS = " << bitset<17>(g_flags_eval_id).to_string() << endl;
      do_vegas_integration(Int_lo_2_6_pdf_x,Int_up_2_6_pdf_x,14,&Integrand_withPDF_2_2,v_state,results[2],errors[2],vp,ip,1,couT);
    }

  ////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 00010  [LO: PHI+INT]
  ////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_B_PHI)
    {
      PS_2_2 ps_gg_tt("gg->tt-bar");
      PS_2_3 ps_t1("decay_ps_t");
      PS_2_3 ps_t2("decay_ps_tb");
      ps_gg_tt.set_child(0,&ps_t1);
      ps_gg_tt.set_child(1,&ps_t2);
      
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(g_flags_eval_v);
      USET_EVAL_ALL(g_flags_eval_id);
      SET_EVAL_B_PHIxQCD(g_flags_eval_v);
      SET_EVAL_B_PHIxPHI(g_flags_eval_v);

      ip.ps     = &ps_gg_tt;
      vp.calls  = g_options.n_calls;
      // integrate LO matrix elements with PDFs
      gsl_monte_vegas_state* v_state = nullptr;

      // only LO histograms get filled
      for (auto e: obs_dist_vector)
	{
	  e.first->SetActive(BOOST_BINARY(00 001 0));
	}

      couT << endl;
      couT << endl << " EVAL_V_FLAGS = " << bitset<17>(g_flags_eval_v).to_string() << endl;
      couT << endl << " EVAL_D_FLAGS = " << bitset<17>(g_flags_eval_id).to_string() << endl;
      do_vegas_integration(Int_lo_2_6_pdf_x,Int_up_2_6_pdf_x,14,&Integrand_withPDF_2_2,v_state,results[2],errors[2],vp,ip,1,couT);
    }

  ////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 00100  [NLO (V): PHI+INT]
  ////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_V)
    {
      PS_2_2 ps_gg_tt("gg->tt-bar");
      PS_2_3 ps_t1("decay_ps_t");
      PS_2_3 ps_t2("decay_ps_tb");
      ps_gg_tt.set_child(0,&ps_t1);
      ps_gg_tt.set_child(1,&ps_t2);
      
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(g_flags_eval_v);
      USET_EVAL_ALL(g_flags_eval_id);
      USET_EVAL_V_PHI0xQCD1(g_flags_eval_v);
      USET_EVAL_V_PHI1xQCD0(g_flags_eval_v);
      SET_EVAL_V_PHIxPHI(g_flags_eval_v);
      //SET_EVAL_V_NF(g_flags_eval_v);
      
      ip.ps     = &ps_gg_tt;
      vp.calls  = g_options.n_calls;
      // integrate LO matrix elements with PDFs
      gsl_monte_vegas_state* v_state = nullptr;
      
      // only NLO virt histograms get filled
      for (auto e: obs_dist_vector)
	{
	  e.first->SetActive(BOOST_BINARY(00 100 0));
	}

      couT << endl;
      couT << endl << " EVAL_V_FLAGS = " << bitset<17>(g_flags_eval_v).to_string() << endl;
      couT << endl << " EVAL_D_FLAGS = " << bitset<17>(g_flags_eval_id).to_string() << endl;
      do_vegas_integration(Int_lo_2_6_pdf_x,Int_up_2_6_pdf_x,14,&Integrand_withPDF_2_2,v_state,results[2],errors[2],vp,ip,1,couT);
    }

  ////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 01000  [NLO (int. Dip.): PHI+INT]
  ////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_D)
    {
      PS_2_2 ps_gg_tt("gg->tt-bar");
      PS_2_3 ps_t1("decay_ps_t");
      PS_2_3 ps_t2("decay_ps_tb");
      ps_gg_tt.set_child(0,&ps_t1);
      ps_gg_tt.set_child(1,&ps_t2);
      
      // set flags to specify which matrix elements will be evaluated
      USET_EVAL_ALL(g_flags_eval_v);
      USET_EVAL_ALL(g_flags_eval_id); 
      SET_EVAL_D_ALL(g_flags_eval_id);
      // set also flags for born contributions that should be evaluated
      SET_EVAL_B_PHIxQCD(g_flags_eval_id);
      SET_EVAL_B_PHIxPHI(g_flags_eval_id);
		  
      ip.ps     = &ps_gg_tt;
      vp.calls  = g_options.n_calls;
      // integrate LO matrix elements with PDFs
      gsl_monte_vegas_state* v_state = nullptr;
      
      // only NLO dip histograms get filled
      for (auto e: obs_dist_vector)
	{
	  e.first->SetActive(BOOST_BINARY(01 000 0));
	}

      couT << endl;
      couT << endl << " EVAL_V_FLAGS = " << bitset<17>(g_flags_eval_v).to_string() << endl;
      couT << endl << " EVAL_D_FLAGS = " << bitset<17>(g_flags_eval_id).to_string() << endl;
      do_vegas_integration(Int_lo_2_6_pdf_x,Int_up_2_6_pdf_x,14,&Integrand_withPDF_2_2,v_state,results[2],errors[2],vp,ip,1,couT);
    }

  ////////////////////////////////////////////////////////////
  // INTEGRATION BLOCK 10000  [NLO (R - uint. Dip.): PHI+INT]
  ////////////////////////////////////////////////////////////
  if (int_flags & I_FLAGS_R)
    {
      PS_2_3 ps_gg_ttg("gg->tt-bar + g");
      PS_2_3 ps_t1("decay_ps_t");
      PS_2_3 ps_t2("decay_ps_tb");
      ps_gg_ttg.set_child(0,&ps_t1);
      ps_gg_ttg.set_child(1,&ps_t2);
      
      // set flags to specify which matrix elements and dipoles will be evaluated
      USET_EVAL_R_ALL(g_flags_eval_r); 
      USET_EVAL_R_FSR_FSR(g_flags_eval_r);
      USET_EVAL_R_ISR_ISR(g_flags_eval_r);
      USET_EVAL_R_FSR_ISR(g_flags_eval_r);
      USET_EVAL_R_FSR_INT(g_flags_eval_r);
      // SET_EVAL_R_FSR_INT(g_flags_eval_r);
      SET_EVAL_R_PHIxPHI_ISR(g_flags_eval_r);
      USET_EVAL_R_PHIxPHI_FSR(g_flags_eval_r);
      // FSR works with arbitrary spin
      
      ip.ps     = &ps_gg_ttg;
      vp.calls  = g_options.n_calls*10;
      gsl_monte_vegas_state* v_state = nullptr;

      // the upper bound of this integration is given by the hadron c.m.e.
      Int_up_2_7_pdf[3] = ip.s_hadr;
      
      // only NLO real histograms get filled
      for (auto e: obs_dist_vector)
	{
	  e.first->SetActive(BOOST_BINARY(10 000 0));
	}

      couT << endl;
      couT << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r).to_string() << endl;
      do_vegas_integration(Int_lo_2_7_pdf,Int_up_2_7_pdf,16,&Integrand_withPDF_2_3,v_state,results[2],errors[2],vp,ip,1,couT);
    }

  double sum = 0.0;
  double err = 0.0;
  for (int i=0;i<5;++i)
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
  couT << "  Born [QCD]         = " << setprecision(prec) << results[0] << " +- " << setprecision(2) << errors[0] << " pb " << endl; 
  couT << "  Born [PHI^2+int.]  = " << setprecision(prec) << results[1] << " +- " << setprecision(2) << errors[1] << " pb " << endl; 
  couT << "  Virtuell           = " << setprecision(prec) << results[2] << " +- " << setprecision(2) << errors[2] << " pb " << endl; 
  couT << "  Dipole int.        = " << setprecision(prec) << results[3] << " +- " << setprecision(2) << errors[3] << " pb " << endl; 
  couT << "  Reell - Dip. uint. = " << setprecision(prec) << results[4] << " +- " << setprecision(2) << errors[4] << " pb " << endl;
  couT << " --------------------------------------------------------------------------------------------- " << endl;
  couT << "  sum                = " << setprecision(prec) << sum << " +- " << setprecision(2) << err << " pb " << endl; 
  couT << " ============================================================================================= " << endl;
      
  
  if (ip.pdf != nullptr) delete ip.pdf;


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  if (g_options.dist)
    {
      bool PLOT_NLO = (int_flags & I_FLAGS_V) || (int_flags & I_FLAGS_D) || (int_flags & I_FLAGS_R);
      
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
      if (int_flags & I_FLAGS_R ) filename << "_R";
      filename << ".root";
      TFile* fresults = new TFile(filename.str().c_str(),"NEW");
            
      double norm[4] = {1.0/Pi2,//results[0],
			1.0/Pi2,//results[0]+results[1],
			0.0,
			1.0/Pi2};//sum};
      // 0 : LO_QCD
      // 1 : LO_PHI
      // 3 : V
      // 4 : D
      // 5 : R

      CanvasPtr c0 = MakeCanvas("Top/Antitop invariant mass distributions");
      DrawDistribution(c0,
		       MttDistributions,
		       &norm[0],
		       string("M_{t#bar{t}} [GeV]"),
		       string("#sigma^{-1} #frac{d#sigma}{dM_{t#bar{t}}} [pb/GeV]"),
		       true,     // WRITE
		       PLOT_NLO);// NLO


      DoubleCanvasPtr cd0 = MakeDoubleCanvas("Top/Antitop transverse momentum distributions");
      DrawTwoDistributions(cd0,
			   PT1Distributions,
			   PT2Distributions,
			   &norm[0],
			   string("p_{T,t} [GeV]"),string("p_{T,#bar{t}} [GeV]"),
			   string("#sigma^{-1} #frac{d#sigma}{dp_{T,t}} [pb/GeV]"),
			   string("#sigma^{-1} #frac{d#sigma}{dp_{T,#bar{t}}} [pb/GeV]"),
			   true,     // WRITE
			   PLOT_NLO);// NLO

      DoubleCanvasPtr cd1 = MakeDoubleCanvas("Top/Antitop pseudorapidity distributions");
      DrawTwoDistributions(cd1,
			   Y1Distributions,
			   Y2Distributions,
			   &norm[0],
			   string("Y_{t}"),string("Y_{#bar{t}}"),
			   string("#sigma^{-1} #frac{d#sigma}{dY_{t}} [pb]"),
			   string("#sigma^{-1} #frac{d#sigma}{dY_{#bar{t}}} [pb]"),
			   true,     // WRITE
			   PLOT_NLO);// NLO
 

      DoubleCanvasPtr cd2 = MakeDoubleCanvas("Lepton/Antilepton transverse momentum distributions");
      DrawTwoDistributions(cd2,
			   PTL1Distributions,
			   PTL2Distributions,
			   &norm[0],
			   string("p_{T,L} [GeV]"),string("p_{T,#bar{L}} [GeV]"),
			   string("#frac{1}{#sigma} #frac{d#sigma}{dp_{T,L}} [pb/GeV]"),
			   string("#frac{1}{#sigma} #frac{d#sigma}{dp_{T,#bar{L}}} [pb/GeV]"),
			   true);

      CanvasPtr c1 = MakeCanvas("Lepton/Antilepton opening angle in lab frame");
      DrawDistribution(c1,
		       PHILLDistributions,
		       &norm[0],
		       string("cos(#phi)"),
		       string("#frac{1}{#sigma} #frac{d#sigma}{d cos(#phi)} [pb]"),
		       true);
      fresults->Close();
      TheApp.Run();
    }
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////

  return 1;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

