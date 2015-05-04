
// cstd lib header files
#include <iostream>
#include <cmath>
#include <bitset>
#include <random>
#include <signal.h>
#include <vector>
using namespace std;

// local header files

#include "../inc/Global.h"
#include "../inc/Functions_Shared.h"
#include "../inc/Makros.h"
#include "../inc/Flags.h"
#include "../inc/Integrands_pp_HX.h"

// LHAPDF header files
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"


#ifdef __cplusplus
extern "C" {
#endif
  extern void qlinit_();
#ifdef __cplusplus
}
#endif


int& int_flags = g_options.int_flags;

// Liebe Gudrun, liebe Anna,
// ich hoffe ihr hattet schöne Ostertage. Meine Brüder und ich waren jetzt auch seit längerer Zeit mal wieder alle zusammen in Römerberg und haben da ein schönes Wochenende verbracht.
// Meine Tage in Aachen sind jetzt tatsächlich gezählt, ich bin nur noch im April vertraglich an die Uni dort gebunden und werden dann Anfang Mai meinen Umzug nach Weinheim komplettieren. Ich hoffe, dass ich auch bis dahin meine Dissertation fertig habe und abgeben kann und dann nur noch einmal zur Verteidigung nach Aachen muss. In Weinheim habe ich ja bereits eine Wohnung zusammen mit meiner Freundin, die dort in der Nähe ihr Referendariat macht, und ich werde dann auch in der Gegend auch Jobsuche gehen. Da ich von der Uni und der Grundlagenforschung nun genug habe werde ich auch mein Glück in der Privatwirtschaft versuchen.



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
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
  if (THDM_1.NBosons() != 1)
    {
      cout << endl << " This program runs with a single Higgs boson! " << endl;
      exit(1);
    }
  
  // create an outpufile
  ofstream log_file;
  if (g_options.logfile) log_file.open(string("results/integration_pp_HX_results_")+=run_id.str());
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
  // some aliases
  double const& mt      = THDM_1.mt();
  double const& mt2     = THDM_1.mt2();
  double const& mScale  = THDM_1.Scale();
  double const& mScale2 = THDM_1.Scale2();
  HPtr Phi1 = THDM_1.GetBoson(0);
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
  double CME = sqrt(ip.s_hadr)*mScale;
  PRINTS(couT,CME);
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // parameters for VEGAS integration routine
  vegas_par vp;
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
    PRINTS(couT,COLL_CUT);
  }
  ////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  Integrator INT(cout);

  Integral Int_2_1(2,IntLimits::Int_lo_2_1_pdf_x,IntLimits::Int_up_2_1_pdf_x);
  Int_2_1.SetIntegrand(&Integrand_2_1_pdf);

  Integral Int_2_2(3,IntLimits::Int_lo_2_2_pdf,IntLimits::Int_up_2_2_pdf);
  Int_2_2.SetIntegrand(&Integrand_2_2_pdf);
  
  
  // mass values to evaluate in main loop
  vector<double> M_values = {100, 125, 150, 200, 250, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600, 700, 800, 900, 1000};
  //vector<double> M_values = {125, 150, 200, 250, 300};
  //vector<double> M_values = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  // output
  /////////////////////////////////////////////////////////////////////////////////////////
  // write mass values to log-file in C-array format
  log_file << endl << "=====================================================" << endl;
  log_file << endl << "double vals_mH[] = {" << endl;
  cout << endl << " Higgs Masses: " << endl << " [";  
  for (auto m = M_values.cbegin(); m != M_values.cend(); ++m)
    {
      log_file << *m;
      cout << *m;
      if (m != --M_values.cend())
	{
	  log_file << ",";
	  cout << ",";
	}
      log_file << endl;
    }
  log_file << "};\n\n";
  cout << "]\n\n";
  
  if (int_flags & BOOST_BINARY(0 000 1)) log_file << " [LO] ";
  if (int_flags & BOOST_BINARY(1 111 0)) log_file << " [NLO] ";
  if (int_flags & BOOST_BINARY(0 001 0)) log_file << " [gg,V] ";
  if (int_flags & BOOST_BINARY(0 010 0)) log_file << " [gg,R] ";
  if (int_flags & BOOST_BINARY(0 100 0)) log_file << " [qq] ";
  if (int_flags & BOOST_BINARY(1 000 0)) log_file << " [qg] ";
  log_file << "\n MU = " << 0.5*g_options.ren_scale << "* mH";
  log_file << "\n\n\n";

  log_file << "#define XSEC(I)  vals_Xsection(I*2+0)\n"; 
  log_file << "#define XERR(I)  vals_Xsection(I*2+1)\n\n";
  
  // write cross-section values and errors to log-file in C-array format
  log_file << endl << "//// result [pb] , error [pb] ////";
  log_file << endl << "double vals_Xsection[] = {" << endl;
  /////////////////////////////////////////////////////////////////////////////////////////


  
  // second loop. now calculate the cross-section to each mass value
  for (auto m = M_values.cbegin(); m != M_values.cend(); ++m)
    {
      // alias Higgs mass
      double const& MH = *m; 
      // reset Higgs mass ...
      Phi1->SetM(MH/mScale);
      // ... as well as alphas and scales (default mu = mH/2, use )
      double MU = 0.5*MH;
      if (g_options.ren_scale!=0.0)
	{
	  MU *= g_options.ren_scale;
	}
      THDM_1.SetAlphaS(ip.pdf->alphasQ(MU));
      THDM_1.SetMUR(MU);
      THDM_1.SetMUF(MU);
      /////////////////////////////////////////////////////////////////////
      // set b-quark MS-bar mass at scale mu [GeV] (O(AlphaS^1)
      // running mass with 5 active flavours, for 6 flavours change exponent to 4.0/7.0)
      // reference value mb(mb) = 4.213 GeV
      double mb     = 4.213/mScale*pow(ip.pdf->alphasQ(MU)/ip.pdf->alphasQ(4.213),12.0/23.0);
      double mb2    = pow(mb,2);
      /////////////////////////////////////////////////////////////////////
      double sigma_tot = 0.0;
      double err_tot   = 0.0;

      // the Higgs mass squared in units of mScale !!!
      double MH2 = pow(MH/mScale,2);

      // set the scaling factor to rescale NLO contributions obtained with eff. ggH vertex
      Phi1->SetFormFactors(MH2,mt2,mb2);
      c_double F_ggH1_s = (-4.0*MH2)*(Phi1->GetFH(0));
      c_double F_ggH1_p = ( 8.0*MH2)*(Phi1->GetFA(0));
      c_double const& FH_eff_1 = Phi1->GetFH(1);
      c_double const& FA_eff_1 = Phi1->GetFA(1);

      // rescaling factor for radiative corrections
      if (g_options.useK)
	{
	  ip.K = (std::norm(F_ggH1_s)+std::norm(F_ggH1_p))/(16.0*MH2*MH2*(std::norm(FH_eff_1)+4.0*std::norm(FA_eff_1)));
	}
      else ip.K = 1.0;
      ////////////////////////////////////////////////////////////
      // INTEGRATION BLOCK 00011  [LO and NLO V (GG)]
      ////////////////////////////////////////////////////////////
      if ((ip.eval_flags = (int_flags & BOOST_BINARY(0 001 1))))
	{
	  ip.SetPS(new PS_2_1(MH2,"pp->H"));
	  // PS_2_1 needs to be set once only when MH2 changes
	  ((PS_2_1*)ip.ps)->set();
	  // adjust # calls
	  vp.calls  = g_options.n_calls/100;
	  
	  INT.Integrate(Int_2_1,ip,vp);
	  sigma_tot += vp.result;
	  err_tot   += vp.error;
	}

      ////////////////////////////////////////////////////////////
      // INTEGRATION BLOCK 11100  [NLO R (GG,QQ,QG)]
      ////////////////////////////////////////////////////////////
      if ((ip.eval_flags = (int_flags & BOOST_BINARY(1 110 0))))
	{
	  // k1 -> Higgsboson, k2 -> massless gluon
	  ip.SetPS(new PS_2_2(MH2,0.0,"pp->Hx"));
	  // adjust # calls
	  vp.calls  = g_options.n_calls/100;

	  INT.Integrate(Int_2_2,ip,vp);
	  sigma_tot += vp.result;
	  err_tot   += vp.error;
	}      

      int prec = (int)log10(sigma_tot/err_tot)+2;
      cout << endl;
      cout << "=== total cross section ============================" << endl;
      cout << " sigma = " << std::setprecision(prec) << sigma_tot << " +- " << std::setprecision(2) << sqrt(err_tot) << " pb " << endl;
      if (int_flags & BOOST_BINARY(0 000 1)) cout << " [LO] ";
      if (int_flags & BOOST_BINARY(1 111 0)) cout << " [NLO] ";
      if (int_flags & BOOST_BINARY(0 001 0)) cout << " [gg,V] ";
      if (int_flags & BOOST_BINARY(0 010 0)) cout << " [gg,R] ";
      if (int_flags & BOOST_BINARY(0 100 0)) cout << " [qq] ";
      if (int_flags & BOOST_BINARY(1 000 0)) cout << " [qg] ";       
      cout << std::setprecision(prec)  << endl;
      cout << " mH [GeV] = " << MH  << endl;
      cout << " mu [GeV] = " << MU << endl;
      cout << " mb [GeV] = " << mb*mScale << endl;
      cout << " AlphaS   = " << THDM_1.AlphaS() << endl;
      cout << " F_s      = " << F_ggH1_s << endl;
      cout << " F_p      = " << F_ggH1_p << endl;
      cout << " fH, fA   = " << FH_eff_1 << ", " << FA_eff_1 << endl;
      cout << " K        = " << ip.K << endl;
      cout << "====================================================" << endl;

      log_file << sigma_tot << ", " << sqrt(err_tot);
      if (m != --M_values.cend()) log_file << ",";
      log_file << endl;
    }
  log_file << "};" << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////





