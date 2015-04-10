
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
  using RunParameters::mScale;
  using RunParameters::mScale2;
  //using RunParameters::mt2;
  using RunParameters::mb;
  using RunParameters::mb2;
  ////////////////////////////////////////////////////////////////////////////////////////////
  // init ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // parse arguments and initialize input parameters
  parse_arguments(argc,argv,g_options);
  ////////////////////////////////////////////////////////////////////////////////
  
  // create identifier from process id and timestamp
  stringstream run_id;
  run_id << "RUN_" << time(0) << "_" << gettid();

  // create an outpufile
  ofstream log_file;
  if (g_options.logfile) log_file.open(string("results/integration_results_")+=run_id.str());
  // the teestream has problems when the file is not opened, so just dump to /dev/null
  else log_file.open("/dev/null");
  
  // create teestream object to tee output to std::couT and file
  teestream couT(std::cout,log_file);
  PRINTS(couT,run_id.str());
  
  // init external libraries
  qlinit_();
#ifdef WITH_NON_FACT_DIAGRAMS
  ltini();
#endif

  /////////////////////////////////////////////////////////////////////////////////////////////
  // input parameters /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////
  // this structure is handed down to the integrand function
  // contains hadronic c.m.e., pointer to distributions, pointer to PDFs
  integrand_par ip;
  ip.s_hadr        = pow((double)g_options.cme*1000.0,2)/mScale2;
  ip.collect_dist  = false;
  ip.pdf           = LHAPDF::mkPDF("CT10nlo", 0);
  // no distributions needed for pp->HX
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
  { // set up Higgs boson parameters
    using namespace HiggsBosons;
    SetPhi1(g_options.mH,
	    0.0,//g_options.GammaH,
	    1.0,//g_options.At,//1.785137274,
	    0.0,//g_options.Bt,
	    1.0,//g_options.Ab,//-.1677099516,
	    0.0 //g_options.Bb
	    );

    PRINTS(couT,Vh*mScale);
    // PRINTS(couT,M_1*mScale);
    PRINTS(couT,G_1*mScale);
    PRINTS(couT,At_1*Vh);
    PRINTS(couT,Bt_1*Vh);
    PRINTS(couT,Ab_1*Vh);
    PRINTS(couT,Bb_1*Vh);
    PRINTS(couT,FH_eff_1/mScale);
    PRINTS(couT,FA_eff_1/mScale);
  }
  ////////////////////////////////////////////////////////////////////////////////
  { // set up run parameters
    using namespace RunParameters;
    // set scales -> will be done in the loop below
    // one or two Higgs bosons?
    TwoHDM = 0;
    PRINTS(couT,mt*mScale);
    PRINTS(couT,AlphaS);
    PRINTS(couT,MUR*mScale);
    PRINTS(couT,MUF*mScale);
  }
  ////////////////////////////////////////////////////////////////////////////////
  { // set technical cuts on soft/collinear 2->3 phase space regions
    using namespace Cuts;
    COLL_CUT = 1.0-pow(10.0,-g_options.tech_cut);
    SOFT_CUT = pow(10.0,-g_options.tech_cut);
    PRINTS(couT,SOFT_CUT);
    PRINTS(couT,COLL_CUT);
  }
  ////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////

  Integrator INT(cout);
  
  Integral Int_2_1(2,IntLimits::Int_lo_2_1_pdf_x,IntLimits::Int_up_2_1_pdf_x);
  Int_2_1.SetIntegrand(&Integrand_2_1_pdf);
  
  // integral over integrated dipoles 
  Integral Int_2_2(3,IntLimits::Int_lo_2_2_pdf,IntLimits::Int_up_2_2_pdf);
  Int_2_2.SetIntegrand(&Integrand_2_2_pdf);


  
  
  // mass values to evaluate in main loop
  vector<double> M_values = {100, 125, 150, 200, 250, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600, 700, 800, 900, 1000};

  // write mass values to log-file in C-array format
  log_file << endl << "=====================================================" << endl;
  log_file << endl << "double vals_mH[] = {" << endl;
  for (auto m = M_values.cbegin(); m != M_values.cend(); ++m)
    {
      log_file << *m;
      if (m != --M_values.cend()) log_file << ",";
      log_file << endl;
    }
  log_file << "};\n\n";

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
  // second loop. now calculate the cross-section to each mass value
  for (auto m = M_values.cbegin(); m != M_values.cend(); ++m)
    {
      using namespace HiggsBosons;
      
      double const& MH = *m; 
      /////////////////////////////////////////////////////////////////////
      HiggsBosons::SetMPhi1(MH); // have to reset mass in every iteration
      /////////////////////////////////////////////////////////////////////
      double MU = 0.5*MH/mScale*g_options.ren_scale;
      RunParameters::SetAlphaS(ip.pdf->alphasQ(MU*mScale));
      RunParameters::SetMUR(MU);
      RunParameters::SetMUF(MU); // have to reset scales in every iteration
      /////////////////////////////////////////////////////////////////////
      // set b-quark MS-bar mass at scale mu [GeV] (O(AlphaS^1)
      // running mass with 5 active flavours, for 6 flavours change exponent to 4.0/7.0)
      // reference value mb(mb) = 4.213 GeV
      mb     = 4.213/mScale*pow(ip.pdf->alphasQ(MU*mScale)/ip.pdf->alphasQ(4.213),12.0/23.0);
      mb2    = pow(mb,2);
      /////////////////////////////////////////////////////////////////////
      double sigma_tot = 0.0;
      double err_tot   = 0.0;

      // the Higgs mass squared in units of mScale !!!
      double MH2 = pow(MH/mScale,2);

      // set the scaling factor to rescale NLO contributions obtained with eff. ggH vertex
      HiggsBosons::set_F_ggH(MH2);
      if (g_options.dist)
	{
	  ip.K = (std::norm(F_ggH1_s)+std::norm(F_ggH1_p))/(16.0*MH2*MH2*(pow(FH_eff_1,2)+4.0*pow(FA_eff_1,2)));
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
	  vp.calls  = g_options.n_calls/10;
	  
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
	  vp.calls  = g_options.n_calls/10;

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
      cout << " mH [GeV] = " << sqrt(M2_1)*mScale  << endl;
      cout << " mu [GeV] = " << MU*mScale << endl;
      cout << " mb [GeV] = " << mb*mScale << endl;
      cout << " AlphaS   = " << RunParameters::AlphaS << endl;
      cout << " F_s      = " << HiggsBosons::F_ggH1_s << endl;
      cout << " F_p      = " << HiggsBosons::F_ggH1_p << endl;
      cout << " fH, fA   = " << HiggsBosons::FH_eff_1 << ", " << HiggsBosons::FA_eff_1 << endl;
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





