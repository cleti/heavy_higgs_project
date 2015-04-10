

#ifndef GLOBAL_H
#define GLOBAL_H

#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <iostream>
#include <streambuf>

#include "Makros.h"

typedef std::complex<double> c_double;
typedef unsigned long ulong;

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Constants {
  /* misc. constants */
  extern const double Pi;
  extern const double gE;
  extern const double Pi2;
  extern const double Pi3;
  extern const double gE2;
  extern const double TwoPi;
  extern const double FourPi;
  extern const double TwoPi2;
  /* conversion factors */
  extern const double CONV_MeV_fm;
  extern const double CONV_GeV2i_mbarn;
  extern const double CONV_GeV2i_pbarn;
  /* first two coefficients of the expansion of (4 pi)^eps/Gamma(1-eps) */
  extern const double C_eps1;
  extern const double C_eps2;
  /* colour factors */
  extern const double CA;
  extern const double TF;
  extern const double CF;
  extern const double CFCA2;
  extern const double Nf;
  extern const double CA2;
  extern const double TF2;
  extern const double CF2;
  extern const double beta0;
  /* initial gg colour and spin average */
  extern const double PREF_GG;
  extern const double PREF_QQ;
  extern const double PREF_QG;
  /* top mass */
  extern const double MT;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace AmpPrefactors {
  /* born */
  extern double PREF_B_PHIxPHI;
  extern double PREF_B_PHIxQCD;
  extern double PREF_B_QCDxQCD;
  extern double PREF_B_QCDxQCD_CF;
  extern double PREF_B_QCDxQCD_CA;
  extern double PREF_B_QCDxQCD_CFCA;
  /* virtual */
  extern double PREF_V;
  extern double PREF_V_CF;
  extern double PREF_V_CA;
  extern double PREF_V_CFCA2;
  extern double PREF_V_Nf;
  extern double PREF_V_CT;
  extern double PREF_V_PHI;
  extern double PREF_V_PHI_CA;
  extern double PREF_V_PHI_CF;
  /* real */
  extern double PREF_R;
  extern double PREF_R_CF;
  extern double PREF_R_CA;
  extern double PREF_R_CFCA2;
  extern double PREF_R_PHI;
  extern double PREF_R_PHI_CA;
  extern double PREF_R_PHI_CF;
  /* UID */
  extern double PREF_UID_TF;
  extern double PREF_UID_CA;
  extern double PREF_UID_CF;
  
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace HiggsBosons {
  // combined VEV
  extern double Vh;
  // PHI1
  // mass & width
  extern double M_1;
  extern double G_1;
  extern double M2_1;
  extern double G2_1;
  // fermion couplings
  extern double At_1;
  extern double Bt_1;
  extern double Ab_1;
  extern double Bb_1;
  // effective gg-Phi couplings
  extern double FH_eff_1;
  extern double FA_eff_1;
  // full 1-loop form factors
  extern c_double F_ggH1_s;
  extern c_double F_ggH1_p;
  
  // PHI1
  // mass & width
  extern double M_2;
  extern double G_2;
  extern double M2_2;
  extern double G2_2;
  // fermion couplings
  extern double At_2;
  extern double Bt_2;
  extern double Ab_2;
  extern double Bb_2;
  // effective gg-Phi couplings
  extern double FH_eff_2;
  extern double FA_eff_2;
  // full 1-loop form factors
  extern c_double F_ggH2_s;
  extern c_double F_ggH2_p;

  // interference terms
  extern double At_fH_re;		
  extern double At_fA_re;		
  extern double Bt_fH_re;		
  extern double Bt_fA_re;
  
  extern double At_fH_im;
  extern double At_fA_im;	
  extern double Bt_fH_im;		
  extern double Bt_fA_im;
  // phi^2
  extern double At2_fH2_De;	
  extern double At2_fA2_De;	
  extern double Bt2_fH2_De;	
  extern double Bt2_fA2_De;
#ifdef WITH_T_SPIN
  extern double At_Bt_fH2_De;	
  extern double At_Bt_fA2_De;
  extern double At_Bt_fH2_DeIM;	
  extern double At_Bt_fA2_DeIM;		
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace RunParameters {
  // flags
  extern unsigned long g_flags_eval_v;
  extern unsigned long g_flags_eval_id;
  extern unsigned long g_flags_eval_r;
  /* mass scale for normalization */
  extern double mScale;
  extern double mScale2;
  // top and gauge masses
  extern double mt;
  extern double mb;
  extern double mW;
  extern double mZ;
  extern double mt2;
  extern double mb2;
  extern double mW2;
  extern double mZ2;
  // widths
  extern double GammaT;
  extern double GammaW;
  extern double GammaT2;
  extern double GammaW2;
  // couplings
  extern double Alpha0;
  extern double Gw2;
  extern double Gw4;
  extern double AlphaS;
  extern double AlphaS2;
  extern double CW2;
  extern double SW2;
  // scales
  extern double MUR;
  extern double MUR2;
  extern double MUF;
  extern double MUF2;
  // deprecated
  extern double LNMU2;
  // conversion factor depends on mass scale for normalization
  extern double CONV_mt2i_pbarn;
  // use both scalar bosons
  extern unsigned TwoHDM;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Cuts {
  extern double COLL_CUT;
  extern double SOFT_CUT;
  extern double IDIP_X_CUT;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// teebuf and teestream
class teebuf: public std::streambuf
{
 public:
  // Construct a streambuf which tees output to both input
  // streambufs.
 teebuf(std::streambuf * sb1, std::streambuf * sb2) :
  sb1(sb1),
  sb2(sb2)
 { }
 private:
  // This tee buffer has no buffer. So every character "overflows"
  // and can be put directly into the teed buffers.
  virtual int overflow(int c)
  {
    if (c == EOF)
      {
	return !EOF;
      }
    else
      {
	int const r1 = sb1->sputc(c);
	int const r2 = sb2->sputc(c);
	return (r1 == EOF || r2 == EOF) ? EOF : c;
      }
  }
    
  // Sync both teed buffers.
  virtual int sync()
  {
    int const r1 = sb1->pubsync();
    int const r2 = sb2->pubsync();
    return (r1 == 0 && r2 == 0) ? 0 : -1;
  }   
 private:
  std::streambuf * sb1;
  std::streambuf * sb2;
};

class teestream : public std::ostream
{
 public:
  // Construct an ostream which tees output to the supplied
  // ostreams.
 teestream(std::ostream & o1, std::ostream & o2) :
  std::ostream(&tbuf),
  tbuf(o1.rdbuf(), o2.rdbuf())
 { }      
 private:
  teebuf tbuf;
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#endif
