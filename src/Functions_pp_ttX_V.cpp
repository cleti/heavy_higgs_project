



// ignore warnings of unused variables
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma message "Note: '-Wunused-variable' disabled in file " __FILE__


#include "../inc/Functions_pp_ttX_V.h"


#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\



//////////////////////////////////////////////////////////////////////////////////////////
// renormalization factors ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// delta Z_M, mass renormalization factor for a quark (MS-bar, 1-loop, 5-flavour)
static inline double ZM (
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
static inline double Z1Qg (
		    double const& MT2,
		    double const& AlphaS,
		    double const& MUR2)
{
  return ZM(MT2,AlphaS,MUR2);
}
// delta Z_2, massive quark wavefunction renormalization factor (MS-bar, 1-loop, 5-flavour)
static inline double Z2 (
		  double const& MT2,
		  double const& AlphaS,
		  double const& MUR2)
{
  return ZM(MT2,AlphaS,MUR2);
}
// delta (Z_2 * Z_M) = delta Z_2 + delta Z_M at one-loop order (MS-bar, 1-loop, 5-flavour)
static inline double Z2M (
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
#include "../amp/virtual/fullSpin/2RE_PHIxQCD_LO_eps0_g4_b.cpp"

// NLO virtual amplitudes
// PHI0 x QCD1
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_SE_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_4Gluon_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Vertex1_eps0_g4.cpp"
#include "../amp/virtual/fullSpin/2RE_HIGGSxQCD_NLO_V_Vertex2_eps0_g4.cpp"
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
  PREFACTORS(hm); // defines references ap and hp
      
  // these have to be adjusted whenever S changes!!!
  if (EVAL_B_PHI(flags))
    {
      hm.SetHiggsPrefactors(ps.get_s(),EFF);
    }
  
  double B_t = 0.0;
  
  //////////////////////////////////////////////////////////////
  // [PHI0]^2: gg initial state
  //////////////////////////////////////////////////////////////
  if (flags & F_EVAL_B_PHIxPHI) 
    {
      B_t += Eval_B_PHIxPHI_withINT12(ps,hm);
#ifdef DEBUG
      CHECKNAN(B_t);
#endif
    }
  
  //////////////////////////////////////////////////////////////
  // 2RE[PHI0 x QCD0]: gg initial state
  //////////////////////////////////////////////////////////////
  if (flags & F_EVAL_B_PHIxQCD) 
    {
      B_t += Eval_B_2PHIxQCD(ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(B_t);
#endif
    }

  //////////////////////////////////////////////////////////////
  // [QCD0]^2: gg initial state
  //////////////////////////////////////////////////////////////
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
  PREFACTORS(hm); // defines references ap and hp 
#ifdef WITH_T_SPIN
  double ret = Eval_B_PHIxPHI(ps,ap,hp);
  // the additional interference terms are only relevant in the polarized amplitudes
  // when models with more than 1 heavy Higgs boson are considered
  if (hm.NBosons()>1)
    {
      ret += Eval_B_PHIxPHI_IM_INTab(ps,ap,hp);
    }
  return ret;
#else
  return Eval_B_PHIxPHI(ps,ap,hp);
#endif
  
}

//////////////////////////////////////////////////////////////
// [QCD0]^2: qq-bar initial state
//////////////////////////////////////////////////////////////
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
  PREFACTORS(hm); // defines references ap and hp 
      
  // these have to be adjusted whenever S changes!!!
  hm.SetHiggsPrefactors(ps.get_s(),1);
  // set scalar integrals (update invariants sometime before this call by Set_Invariants())
  Set_SI(ps.get_msq(0), // assuming that get_msq(1) is the same value
	 ps.get_rs(),
	 ps.get_t11(),
	 ps.get_t12(),
	 hm.MUR2(),
	 flags);

  double res = 0.0;

  
  //////////////////////////////////////////////////////////////
  // 2RE[PHI0 x QCD1] interference terms
  // labels according to the conventions introduced in my thesis 
  //////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (flags & F_EVAL_V_PHI0xQCD1)
    {
      if (flags & F_EVAL_V_SE)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // self-energy insertion diagram (+ren. counter term)
	  res += Eval_V_SE (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}
#endif
      
      
#ifdef DEBUG
      if (flags & F_EVAL_V_4G)
	{
#endif         
	  //////////////////////////////////////////////////////////
	  // 4g vertex diagram
	  res += Eval_V_4G (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
      	}
#endif
      
      
#ifdef DEBUG
      if (flags & F_EVAL_V_V1)
	{
#endif
#ifdef WITH_T_SPIN
	  //////////////////////////////////////////////////////////
	  // Qg vertex correction (+ren. counter term)
	  // full top spin dependency: upper vertex != lower vertex
	  res += Eval_V_V1 (ps,ap,hp);
	  res += Eval_V_V2 (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#else
	  //////////////////////////////////////////////////////////
	  // Qg vertex correction (+ren. counter term)
	  // summed over top spin: upper vertex = lower vertex
	  res += 2.0 * Eval_V_V1 (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#endif
#ifdef DEBUG
	  CHECKNAN(res);
      	}
#endif
      
     
#ifdef DEBUG
      if (flags & F_EVAL_V_B1)
	{
#endif         
	  //////////////////////////////////////////////////////////
	  // box diagram 1
	  res += Eval_V_B1 (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
      	}
#endif
      
     
#ifdef DEBUG
      if (flags & F_EVAL_V_B2)
	{
#endif         
	  //////////////////////////////////////////////////////////
	  // box diagram 2
	  res += Eval_V_B2 (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
      	}
#endif
      
      
#ifdef DEBUG
      if (flags & F_EVAL_V_B3)
	{
#endif         
	  //////////////////////////////////////////////////////////
	  // box diagram 3
	  res += Eval_V_B3 (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
      	}
    }
#endif


  //////////////////////////////////////////////////////////////
  // 2RE[PHI1 x QCD0] interference terms
  //////////////////////////////////////////////////////////////
#ifdef DEBUG  
  if (flags & F_EVAL_V_PHI1xQCD0)
    {
#endif
      res += Eval_V_2RE_PHI1xQCD0(ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif


  //////////////////////////////////////////////////////////////
  // 2RE[PHI1 x PHI0] 
  //////////////////////////////////////////////////////////////  
#ifdef DEBUG 
  if (flags & F_EVAL_V_PHIxPHI)
    {
#endif      
      res += Eval_V_PHIxPHI(ps,ap,hp);
#ifdef WITH_T_SPIN
      if (hm.NBosons()>1)
	{
	  res += Eval_V_PHIxPHI_IM_INTab(ps,ap,hp);
	}
#endif
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif  
  return res;
}





