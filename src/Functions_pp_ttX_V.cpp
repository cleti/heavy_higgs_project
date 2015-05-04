
#include "../inc/Functions_pp_ttX_V.h"


#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\



//////////////////////////////////////////////////////////////////////////////////////////
// renormalization factors ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// delta Z_M, mass renormalization factor for a quark (MS-bar, 1-loop, 5-flavour)
static double ZM (
		  double const& MT2,
		  double const& AlphaS,
		  double const& MUR2)
{
  using namespace Constants;
  static const int I = 0;
  static const double ZERO = 0.0;
  static double res = 0.0;
  static double MT2_t = 0.0;
  static double LNMUR2 = 0.0;
  
  if (MT2_t != MT2)
    {      
      LNMUR2 = (RE(qli2_(&MT2,&ZERO,&MT2,&MUR2,&I))-2.0);//  (RE(I2_MT2_0_MT2_MU2_0)-2.0) = log(MUR2/MT2)
      MT2_t = MT2;
      res = (AlphaS/FourPi*CF)*(-3.0*LNMUR2-4.0);
    } 
  return res;
}
// delta Z_1, massive quark gluon vertex renormalization factor (MS-bar, 1-loop, 5-flavour)
inline double Z1Qg (
		    double const& MT2,
		    double const& AlphaS,
		    double const& MUR2)
{
  return ZM(MT2,AlphaS,MUR2);
}
// delta Z_2, massive quark wavefunction renormalization factor (MS-bar, 1-loop, 5-flavour)
inline double Z2 (
		  double const& MT2,
		  double const& AlphaS,
		  double const& MUR2)
{
  return ZM(MT2,AlphaS,MUR2);
}
// delta (Z_2 * Z_M) = delta Z_2 + delta Z_M at one-loop order (MS-bar, 1-loop, 5-flavour)
inline double Z2M (
		   double const& MT2,
		   double const& AlphaS,
		   double const& MUR2)
{
  return Z2(MT2,AlphaS,MUR2)+ZM(MT2,AlphaS,MUR2);
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
///////////////////////////////////////////////////////////////////////////////////////////////
// LO amplitudes
// PHI^2
#include "../amp/virtual/fullSpin/PHIxPHI_LO_eps0_g4.cpp"
//#include "../amp/virtual/fullSpin/PHIxPHI_LO_INT12_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/PHIxPHI_LO_IM_INTab_eps0_g4.cpp"
// QCD^2
#include "../amp/virtual/fullSpin/QCDxQCD_LO_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/QCDxQCD_LO_eps0_g4_qq.cpp"
// PHI x QCD
#include "../amp/virtual/fullSpin/2RE_PHIxQCD_LO_eps0_g4.cpp"

// NLO virtual amplitudes
// PHI0 x QCD1
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_SE_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_4Gluon_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Vertex1_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Box1_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Box2_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Box3_eps0_g4.cpp"

// #ifdef WITH_NON_FACT_DIAGRAMS
// #include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_NF_1_eps0_g4.cpp"
// #include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_NF_2_eps0_g4.cpp"
// #include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_NF_3_eps0_g4.cpp"
// #endif

// // PHI1 x QCD0
#include "../amp/virtual/fullSpin/2RE_PHI1xQCD0_NLO_V_eps0_g4.cpp"
// // PHI^2
#include "../amp/virtual/fullSpin/PHIxPHI_NLO_V_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/PHIxPHI_NLO_V_IM_INTab_eps0_g4.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////////
// LO amplitudes
// PHI^2
#include "../amp/virtual/noSpin/PHIxPHI_LO_eps0_g4.cpp"
// QCD^2
#include "../amp/virtual/noSpin/QCDxQCD_LO_eps0_g4_gg.cpp"
#include "../amp/virtual/noSpin/QCDxQCD_LO_eps0_g4_qq.cpp"
// PHI x QCD
#include "../amp/virtual/noSpin/2RE_PHIxQCD_LO_eps0_g4.cpp"

// NLO virtual amplitudes
// PHI0 x QCD1
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_SE_eps0_g4.cpp"
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_4Gluon_eps0_g4.cpp"
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_Vertex1_eps0_g4.cpp"
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_Box1_eps0_g4.cpp"
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_Box2_eps0_g4.cpp"
#include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_Box3_eps0_g4.cpp"

// #ifdef WITH_NON_FACT_DIAGRAMS
// #include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_NF_1_eps0_g4.cpp"
// #include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_NF_2_eps0_g4.cpp"
// #include "../amp/virtual/noSpin/2RE_HIGGSxQCD_NLO_V_NF_3_eps0_g4.cpp"
// #endif

// // PHI1 x QCD0
#include "../amp/virtual/noSpin/2RE_PHI1xQCD0_NLO_V_eps0_g4.cpp"
// // PHI^2
#include "../amp/virtual/noSpin/PHIxPHI_NLO_V_eps0_g4.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////////////////////

// BORN
// squared QCD gg->tt amplitude +
// squared gg->PHI->tt amplitude +
// 2RE[QCDxPHI] interference terms
double Eval_B(
	      const PS_2_2& ps,
	      HiggsModel& hm,
	      const ulong& flags,
	      unsigned EFF)
{
  PREFACTORS(hm); // defines ap and hp
      
  // these have to be adjusted whenever S changes!!!
  if (EVAL_B_PHI(flags))
    {
      hm.SetHiggsPrefactors(ps.get_s(),EFF);
      //HiggsBosons::ResetHiggsPrefactors(ps.get_s(),EFF);
      // PRINT(hp.At2_fH2_De);
      // PRINT(hp.Bt2_fH2_De);
      // PRINT(hp.At2_fA2_De);
      // PRINT(hp.Bt2_fA2_De);
      // PRINT(hp.At_Bt_fH2_De);
      // PRINT(hp.At_Bt_fA2_De);
      // exit(1);
    }

       
  double B_t = 0.0;
  if (flags & F_EVAL_B_PHIxPHI) 
    {
      B_t += Eval_B_PHIxPHI_withINT12(ps,hm);
#ifdef DEBUG
      CHECKNAN(B_t);
#endif
    }
  if (flags & F_EVAL_B_PHIxQCD) 
    {
      B_t += Eval_B_2PHIxQCD(ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(B_t);
#endif
    }
  if (flags & F_EVAL_B_QCDxQCD) 
    {
      B_t += Eval_B_QCDxQCD_GG(ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(B_t);
#endif
    }
  return B_t;
}

double Eval_B_PHIxPHI_withINT12(
				const PS_2_2& ps,
				const HiggsModel& hm)
{
  PREFACTORS(hm); // defines ap and hp 
  double ret = Eval_B_PHIxPHI(ps,ap,hp);
#ifdef WITH_T_SPIN
  // the additional interference terms are only relevant in the polarized amplitudes
  // when models with more than 1 heavy Higgs boson are considered
  if (hm.NBosons()>1)
    {
      ret -= Eval_B_PHIxPHI_IM_INTab(ps,ap,hp);
    }
#endif  
  return ret;
}

// squared QCD qq->tt amplitude
double Eval_B_QQ(
		 const PS_2_2& ps,
		 HiggsModel& hm)
{
  PREFACTORS(hm); // defines ap and hp
  return Eval_B_QCDxQCD_QQ(ps,ap,hp);
}


// VIRTUAL CORRECTIONS (finite part)
double Eval_V(
	      const PS_2_2& ps,
	      HiggsModel& hm,
	      const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp 
      
  // these have to be adjusted whenever S changes!!!
  hm.SetHiggsPrefactors(ps.get_s(),1);
  // set scalar integrals (update invariants sometime before this call by Set_Invariants())
  Set_SI(ps,hm.MUR2(),flags);

  double res = 0.0;
  // add contributions of individual diagrams according to flags in flags
  if (flags & F_EVAL_V_PHI0xQCD1)
    {
      // 2RE[PHI0 x QCD1] interference terms classified by QCD1 diagram type
      // self-energy diagram including ren. counter term
      res += Eval_V_SE (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
      // 4g vertex diagram
      res += Eval_V_4G (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
      // qg veretx correction (upper and lower vertex) including ren. counter term
      res += 2.0*Eval_V_V1 (ps,ap,hp); 
#ifdef DEBUG
      CHECKNAN(res);
#endif
      // box 1
      res += Eval_V_B1 (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
      // box 2
      res += Eval_V_B2 (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
      // box 3
      res += Eval_V_B3 (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  // 2RE[PHI1 x QCD0] interference terms
  if (flags & F_EVAL_V_PHI1xQCD0)
    {
      res += Eval_V_2RE_PHI1xQCD0(ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  // 2RE[PHI1xPHI0] interference terms
  if (flags & F_EVAL_V_PHIxPHI)
    {
      res += Eval_V_PHIxPHI(ps,ap,hp);
#ifdef WITH_T_SPIN
      if (hm.NBosons()>1)
	{
	  res -= Eval_V_PHIxPHI_IM_INTab(ps,ap,hp);
	}
#endif
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
    return res;
}





