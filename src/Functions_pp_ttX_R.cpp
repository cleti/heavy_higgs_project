
#include "../inc/Functions_pp_ttX_R.h"



#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\

///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
///////////////////////////////////////////////////////////////////////////////////////////////////
// PHI^2 amplitudes
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_qq.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_qg.cpp"
// // test: my own PHI^2 amplitudes 
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_v2_gg.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_v2_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_INT12_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_INT12_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_IM_INTab_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_IM_INTab_gg.cpp"
// interference terms 2Re[QCDxPHI]
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIINT_gg.cpp"
// // #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIINT_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIINT_gg.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////////////
// PHI^2 amplitudes
// [gg->PHIg->ttg]
#include "../amp/real/noSpin/PHIxPHI_NLO_R_ISR_gg.cpp"
#include "../amp/real/noSpin/PHIxPHI_NLO_R_FSR_gg.cpp"
// [qq->PHIq->ttq]
#include "../amp/real/noSpin/PHIxPHI_NLO_R_qq.cpp"
// [qg->PHIg->ttg]
#include "../amp/real/noSpin/PHIxPHI_NLO_R_qg.cpp"

// interference terms 2Re[QCDxPHI]
// [gg->ttg]
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIISR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIFSR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIINT_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIISR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIFSR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIINT_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIISR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIFSR_gg.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIINT_gg.cpp"
// [qg->ttq]
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIINT_QG.cpp"
// [qq->ttg]
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIINT_QQ.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////




// Real corrections: GG initial state
double Eval_R_GG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp

  double res = 0.0;
  //////////////////////////////////////////////////////////////////////////
  // call before all PHI-FSR contributions with p1+p2 in Phi-denominator
  //////////////////////////////////////////////////////////////////////////
  double s12 = ps.get_s();
  hm.SetHiggsPrefactors(s12,1);
  // interference terms QCDxPHI
#ifndef WITH_T_SPIN
  if (flags & F_EVAL_R_FSR_FSR) // DIV !
    {
      res += Eval_R_FSR_FSR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_ISR_FSR) // DIV !
    {
      res += Eval_R_ISR_FSR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_INT_FSR)
    {
      res += Eval_R_INT_FSR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
#endif // WITH_T_SPIN
  // PHI^2 terms
  if (flags & F_EVAL_R_PHIxPHI_FSR)
    {
      res += Eval_R_PHIxPHI_FSR(ps,ap,hp);

#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (RunParameters::TwoHDM==1)
	{ 
	  res -= Eval_R_PHIxPHI_FSR_IM_INTab(ps,ap,hp);
	}
#endif

#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  //////////////////////////////////////////////////////////////////////////
  // call before all PHI-ISR / INT contributions with k1+k2 in Phi-denominator
  //////////////////////////////////////////////////////////////////////////
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  // interference terms QCDxPHI
    
#ifndef WITH_T_SPIN
  if (flags & F_EVAL_R_ISR_ISR) // DIV !
    {
      res += Eval_R_ISR_ISR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_INT_ISR)
    {
      res += Eval_R_INT_ISR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_FSR_ISR) // DIV ! (corresponds to non-factorizable virtual amplitudes)
    {
      res += Eval_R_FSR_ISR (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_ISR_INT)
    {
      res += Eval_R_ISR_INT (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_FSR_INT) // (corresponds to non-factorizable virtual amplitudes)
    {
      res += Eval_R_FSR_INT (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  if (flags & F_EVAL_R_INT_INT)
    {
      res += Eval_R_INT_INT (ps,ap,hp);
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
#endif // WITH_T_SPIN
  //Phi^2 terms
  if (flags & F_EVAL_R_PHIxPHI_ISR)
    {
      res += Eval_R_PHIxPHI_ISR(ps,ap,hp);

#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (RunParameters::TwoHDM==1)
	{
	  res -= Eval_R_PHIxPHI_ISR_IM_INTab(ps,ap,hp);
	}
#endif

#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  return res;
}




// Real corrections: QQ initial state
double Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp

  double res = 0.0;

  //////////////////////////////////////////////////////////////////////////
  // call before all PHI-ISR / INT contributions with k1+k2 in Phi-denominator
  //////////////////////////////////////////////////////////////////////////
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  // interference terms QCDxPHI
  if (flags & F_EVAL_R_PHIxQCD_QQ)
    {
      res += Eval_R_FSR_INT_QQ(ps,ap,hp); 
    }
  //Phi^2 terms
  if (flags & F_EVAL_R_PHIxPHI_QQ)
    {
      res += Eval_R_PHIxPHI_QQ(ps,ap,hp);
      
#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (RunParameters::TwoHDM==1)
	{
	  // res -= Eval_R_PHIxPHI_QQ_IM_INTab(ps);
	}
#endif
      
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  return res;
}






// Real corrections: QG initial state
double Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp

  double res = 0.0;

  //////////////////////////////////////////////////////////////////////////
  // call before all PHI-ISR / INT contributions with k1+k2 in Phi-denominator
  //////////////////////////////////////////////////////////////////////////
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  // interference terms QCDxPHI
  if (flags & F_EVAL_R_PHIxQCD_QG)
    {
      res += Eval_R_INT_INT_QG(ps,ap,hp);
    }
  //Phi^2 terms
  if (flags & F_EVAL_R_PHIxPHI_QG)
    {
      res += Eval_R_PHIxPHI_QG(ps,ap,hp);
      
#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (RunParameters::TwoHDM==1)
	{
	  // res -= Eval_R_PHIxPHI_QG_IM_INTab(ps);
	}
#endif
      
#ifdef DEBUG
      CHECKNAN(res);
#endif
    }
  return res;
}





