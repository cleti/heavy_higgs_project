

// ignore warnings of unused variables
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma message "Note: '-Wunused-variable' disabled in file " __FILE__



#include "../inc/Functions_pp_ttX_R.h"



#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	

///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
///////////////////////////////////////////////////////////////////////////////////////////////////
// PHI^2 amplitudes
// [gg -> PHI g-> tt g]
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_gg.cpp"
// [qq -> PHI q-> tt q]
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_qq.cpp"
// [qg -> PHI g-> tt g]
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_qg.cpp"
// PHI1 x PHI2 interference terms
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_IM_INTab_gg.cpp"
#include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_IM_INTab_gg.cpp"
// // test: my own PHI^2 amplitudes 
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_v2_gg.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_v2_gg.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_ISR_INT12_gg.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_NLO_R_FSR_INT12_gg.cpp"
// interference terms 2Re[QCDxPHI]
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDISR_PHIINT_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIINT_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIISR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIFSR_gg.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIINT_gg.cpp"
// [qg->ttq]
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDINT_PHIINT_QG.cpp"
// [qq->ttg]
// #include "../amp/real/fullSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIINT_QQ.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////////////
// PHI^2 amplitudes
// [gg -> PHI g-> tt g]
#include "../amp/real/noSpin/PHIxPHI_NLO_R_ISR_gg.cpp"
#include "../amp/real/noSpin/PHIxPHI_NLO_R_FSR_gg.cpp"
// [qq -> PHI q-> tt q]
#include "../amp/real/noSpin/PHIxPHI_NLO_R_qq.cpp"
// [qg -> PHI g-> tt g]
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

  //////////////////////////////////////////////////////////////
  // 2RE[PHI_R x QCD_R] interference terms
  // no initial/internal gluon radiation in the PHI diagram
  //  -> the PHI propagator is evaluated at s = (p1+p2)^2
  // labels according to the conventions introduced in my thesis
  double const& s12 = ps.get_s();
  hm.SetHiggsPrefactors(s12,1);
  //////////////////////////////////////////////////////////////

#ifndef WITH_T_SPIN
#ifdef DEBUG
  if (flags & F_EVAL_R_GG) 
    {
      
      if (flags & F_EVAL_R_FSR_FSR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // final-final (divergent)
	  res += Eval_R_FSR_FSR (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}
      
      if (flags & F_EVAL_R_ISR_FSR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // initial-final (divergent)
	  res += Eval_R_ISR_FSR (ps,ap,hp); 
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_INT_FSR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // internal-final
	  res += Eval_R_INT_FSR (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}
      
    }
#endif // DEBUG
#endif // WITH_T_SPIN



#ifdef DEBUG
  // PHI^2 terms
  if (flags & F_EVAL_R_PHIxPHI_FSR)
    {
#endif
      //////////////////////////////////////////////////////////
      // [PHI_R]^2 final state radiation
      //  -> the PHI propagator is evaluated at s = (p1+p2)^2
      res += Eval_R_PHIxPHI_FSR(ps,ap,hp);
      //////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (hm.NBosons()>1)
	{
	  res += Eval_R_PHIxPHI_FSR_IM_INTab(ps,ap,hp);
	}
#endif
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif



  //////////////////////////////////////////////////////////////
  // 2RE[QCD_R x PHI_R] interference terms
  // with initial/internal gluon radiation in the PHI diagram
  //  -> the PHI propagator is evaluated at s = (k1+k2)^2
  // labels according to the conventions introduced in my thesis
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  //////////////////////////////////////////////////////////////
#ifndef WITH_T_SPIN
#ifdef DEBUG
  if (flags & F_EVAL_R_GG) 
    {
      
      if (flags & F_EVAL_R_ISR_ISR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // initial-initial (divergent)
	  res += Eval_R_ISR_ISR (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_INT_ISR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // internal-initial
	  res += Eval_R_INT_ISR (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_FSR_ISR)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // final-initial (divergent)
	  // corresponds to non-factorizable virtual amplitudes
	  res += Eval_R_FSR_ISR (ps,ap,hp); 
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_ISR_INT)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // initial-internal
	  res += Eval_R_ISR_INT (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_FSR_INT)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // final-internal
	  // corresponds to non-factorizable virtual amplitudes
	  res += Eval_R_FSR_INT (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}

      if (flags & F_EVAL_R_INT_INT)
	{
#endif
	  //////////////////////////////////////////////////////////
	  // internal-internal
	  res += Eval_R_INT_INT (ps,ap,hp);
	  //////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(res);
	}
      
    }
#endif // DEBUG
#endif // WITH_T_SPIN


  
#ifdef DEBUG
  if (flags & F_EVAL_R_PHIxPHI_ISR)
    {
#endif
      //////////////////////////////////////////////////////////
      // [PHI_R]^2 initial state and internal radiation
      //  -> the PHI propagator is evaluated at s = (k1+k2)^2
      res += Eval_R_PHIxPHI_ISR(ps,ap,hp);
      //////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
      // these additional interference terms are only needed when spins are involved
      if (hm.NBosons()>1)
	{
	  res += Eval_R_PHIxPHI_ISR_IM_INTab(ps,ap,hp);
	}
#endif
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif
  
  return res;
}





double _Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp

  double res = 0.0;

#ifndef WITH_T_SPIN
#ifdef DEBUG
  if (flags & F_EVAL_R_PHIxQCD_QQ)
    {
#endif
      //////////////////////////////////////////////////////////////
      // 2RE[QCD_R x PHI_R]   
      res += Eval_R_FSR_INT_QQ(ps,ap,hp);
      //////////////////////////////////////////////////////////////
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif // DEBUG
#endif // WITH_T_SPIN
  
#ifdef DEBUG
  if (flags & F_EVAL_R_PHIxPHI_QQ)
    {
#endif
      //////////////////////////////////////////////////////////////
      // [PHI_R]^2 
      res += Eval_R_PHIxPHI_QQ(ps,ap,hp);
      //////////////////////////////////////////////////////////////  

      // why is this commented out???
      // #ifdef WITH_T_SPIN
      //       if (hm.NBosons()>1)
      // 	{
      // 	  res += Eval_R_PHIxPHI_QQ_IM_INTab(ps);
      // 	}
      // #endif
      
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif
  return res;
}


double Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  //////////////////////////////////////////////////////////////
  // real corrections with qq-bar initial state
  //  -> the PHI propagator is evaluated at s = (k1+k2)^2
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  //////////////////////////////////////////////////////////////
  return _Eval_R_QQ(ps,hm,flags);
}





double _Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp

  double res = 0.0;

#ifndef WITH_T_SPIN
#ifdef DEBUG
  if (flags & F_EVAL_R_PHIxQCD_QG)
    {
#endif
      //////////////////////////////////////////////////////////////
      // 2RE[QCD_R x PHI_R] (divergent)
      res += Eval_R_INT_INT_QG(ps,ap,hp);
      //////////////////////////////////////////////////////////////  
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif  // DEBUG
#endif  // WITH_T_SPIN

  
#ifdef DEBUG
  if (flags & F_EVAL_R_PHIxPHI_QG)
    {
#endif
      //////////////////////////////////////////////////////////////
      // [PHI_R]^2 (divergent)
      res += Eval_R_PHIxPHI_QG(ps,ap,hp);
      ////////////////////////////////////////////////////////////// 

      // why is this commented out???
      // #ifdef WITH_T_SPIN
      //       // these additional interference terms are only needed when spins are involved
      //       if (hm.NBosons()>1)
      // 	{
      // 	  // res -= Eval_R_PHIxPHI_QG_IM_INTab(ps);
      // 	}
      // #endif
      
#ifdef DEBUG
      CHECKNAN(res);
    }
#endif
  return res;
}



double Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags)
{
  //////////////////////////////////////////////////////////////
  // real corrections with qg initial state
  //  -> the PHI propagator is evaluated at s = (k1+k2)^2
  double S12 = 2.0*(hm.mt2() + sp(ps.k1(),ps.k2()));
  hm.SetHiggsPrefactors(S12,1);
  //////////////////////////////////////////////////////////////
  return _Eval_R_QG(ps,hm,flags);
}

