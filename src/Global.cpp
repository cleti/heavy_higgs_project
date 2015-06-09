

#include "../inc/Global.h"

  



#include <cmath>
using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Constants {
  /* misc. constants */
  const double Pi	   = 3.141592653589793;//M_PI;
  const double gE          = 0.5772156649015329;
  const double TwoPi	   = 2.0*Pi;
  const double FourPi	   = 4.0*Pi;
  const double Pi2         = pow(Pi,2);
  const double Pi3         = pow(Pi,3);
  const double gE2         = pow(gE,2);
  const double TwoPi2      = pow(TwoPi,2);
  /* conversion factors */
  const double CONV_MeV_fm      = 197.3269718;
  const double CONV_GeV2i_mbarn = 0.389379338;
  const double CONV_GeV2i_pbarn  = CONV_GeV2i_mbarn*std::pow(10.0,9);
  /* first two coefficients of the expansion of (4 pi)^eps/Gamma(1-eps) */
  const double C_eps1 = 0.0;//log(FourPi)-gE;
  const double C_eps2 = 0.0;//0.5*pow(log(FourPi),2)-Pi2/12.0+0.5*gE2-gE*log(FourPi);
  /* colour factors */
  const double CA	 = 3.0;
  const double TF	 = 0.5;
  const double CF	 = (pow(CA,0.2e1)-0.1e1)/(0.2e1*CA);
  const double CFCA2	 = CF-0.5*CA;
  const double Nf	 = 5.0; // # light quarks in loops -> running coupling
  const double CA2	 = pow(CA,2);
  const double TF2	 = pow(TF,2);
  const double CF2	 = pow(CF,2);
  const double beta0     = (11./12.)*CA-Nf/6.0;
  /* initial gg colour and spin average */
  const double PREF_GG = 4.0/(4.0*pow(2.0*CF*CA,2));
  const double PREF_QQ = 4.0/(4.0*CA2);
  const double PREF_QG = 4.0/(4.0*2.0*CF*CA2);

  const double BR_TT_LL = 24./81.;

  const double kappa_p =  0.9985;
  const double kappa_m = -0.9985;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Cuts {
  double COLL_CUT   = 1.0-1.0e-7;
  double SOFT_CUT   = 1e-7;
  double IDIP_X_CUT = 1.0-1.0e-7;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
