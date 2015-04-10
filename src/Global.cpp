

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
  const double PREF_GG = 1.0/(4.0*pow(2.0*CF*CA,2));
  const double PREF_QQ = 1.0/(4.0*CA2);
  const double PREF_QG = 1.0/(4.0*2.0*CF*CA2);
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace AmpPrefactors {
  /* born */
  double PREF_B_PHIxPHI = 0.0;
  double PREF_B_PHIxQCD = 0.0;
  double PREF_B_QCDxQCD = 0.0;
  double PREF_B_QCDxQCD_CF = 0.0;
  double PREF_B_QCDxQCD_CA = 0.0;
  double PREF_B_QCDxQCD_CFCA = 0.0;
  /* virtual */
  double PREF_V    = 0.0;
  double PREF_V_CF = 0.0;
  double PREF_V_CA = 0.0;
  double PREF_V_CFCA2 = 0.0;
  double PREF_V_Nf = 0.0;
  double PREF_V_CT = 0.0;
  double PREF_V_PHI = 0.0;
  double PREF_V_PHI_CA = 0.0;
  double PREF_V_PHI_CF = 0.0;
  /* real */
  double PREF_R = 0.0;
  double PREF_R_CF = 0.0;
  double PREF_R_CA = 0.0;
  double PREF_R_CFCA2 = 0.0;
  double PREF_R_PHI = 0.0;
  double PREF_R_PHI_CA = 0.0;
  double PREF_R_PHI_CF = 0.0;
  /* UID */
  double PREF_UID_TF = 0.0;
  double PREF_UID_CA = 0.0;
  double PREF_UID_CF = 0.0;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace HiggsBosons {
  using namespace RunParameters;
  // combined VEV
  double Vh     = 246.0/mScale;
  // PHI1
  // mass & width
  double M_1  = 500.0/mScale;
  double G_1  = 3.245439053359254e+01/mScale;
  double M2_1 = pow(M_1,2);
  double G2_1 = pow(G_1,2);
  // fermion couplings
  double At_1 = 1.0/Vh;
  double Bt_1 = 1.0/Vh;
  double Ab_1 = 1.0/Vh;
  double Bb_1 = 1.0/Vh;
  // effective gg-Phi couplings
  double FH_eff_1 =  At_1/12.0;
  double FA_eff_1 = -Bt_1/16.0;
  // full 1-loop form factors
  c_double F_ggH1_s = 0.0;
  c_double F_ggH1_p = 0.0;

  // PHI1
  // mass & width
  double M_2  = 550.0/mScale;
  double G_2  = 25.0/mScale;
  double M2_2 = pow(M_2,2);
  double G2_2 = pow(G_2,2);
  // fermion couplings
  double At_2 = 1.0/Vh;
  double Bt_2 = 1.0/Vh;
  double Ab_2 = 1.0/Vh;
  double Bb_2 = 1.0/Vh;
  // effective gg-Phi couplings
  double FH_eff_2 =  At_2/12.0;
  double FA_eff_2 = -Bt_2/16.0;
  // full 1-loop form factors
  c_double F_ggH2_s = 0.0;
  c_double F_ggH2_p = 0.0;
  
  // interference terms
  double At_fH_re = 0.0;		
  double At_fA_re = 0.0;		
  double Bt_fH_re = 0.0;		
  double Bt_fA_re = 0.0;
  
  double At_fH_im = 0.0;
  double At_fA_im = 0.0;	
  double Bt_fH_im = 0.0;		
  double Bt_fA_im = 0.0;
  // phi^2
  double At2_fH2_De = 0.0;	
  double At2_fA2_De = 0.0;	
  double Bt2_fH2_De = 0.0;	
  double Bt2_fA2_De = 0.0;
#ifdef WITH_T_SPIN
  double At_Bt_fH2_De = 0.0;	
  double At_Bt_fA2_De = 0.0;
  double At_Bt_fH2_DeIM = 0.0;	
  double At_Bt_fA2_DeIM = 0.0;
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace RunParameters {
  // flags
  unsigned long g_flags_eval_v  = 0;
  unsigned long g_flags_eval_id = 0;
  unsigned long g_flags_eval_r  = 0;
  /* mass scale for normalization */
  double mScale  = 173.5;
  double mScale2 = pow(mScale,2);
  // top and gauge masses
  double mt = 173.5/mScale;
  double mb = 4.7/mScale;
  double mW = 80.385/mScale;
  double mZ = 91.19/mScale;
  double mt2 = pow(mt,2);
  double mb2 = pow(mb,2);
  double mW2 = pow(mW,2);
  double mZ2 = pow(mZ,2);
  // widths
  double GammaT = 1.35/mScale;
  double GammaW = 2.08/mScale;
  double GammaT2 = pow(GammaT,2);
  double GammaW2 = pow(GammaW,2);
  // couplings
  double Alpha0  = 1.0/132.351;
  double CW2     = mW2/mZ2;
  double SW2     = 1.0-CW2;
  double Gw2     = Constants::FourPi*Alpha0/(2.0*SW2);
  double Gw4     = pow(Gw2,2);
  double AlphaS  = 0.1184;
  double AlphaS2 = pow(AlphaS,2);
  // scales
  double MUR  = 1.0;
  double MUR2 = 1.0;
  double MUF  = 1.0;
  double MUF2 = 1.0;
  // deprecated
  double LNMU2= log(MUR2);
  // conversion factor depends on mass scale for normalization
  double CONV_mt2i_pbarn = Constants::CONV_GeV2i_pbarn/mScale2;
  // 0: only one neutral Hig
  // 1: include Phi1-Phi2 interferences via adjustment of Prefactors
  // 2: include Phi1-Phi2 interferences explicitly
  unsigned TwoHDM = 0;
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
