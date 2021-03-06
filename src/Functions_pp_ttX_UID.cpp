


// ignore warnings of unused variables
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma message "Note: '-Wunused-variable' disabled in file " __FILE__



#include "../inc/Functions_pp_ttX_UID.h"




#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\


// c.f. (5.42) in hep-ph/0201036v1
static inline double X (const FV&p_i, const FV&p_j, const FV&p_a)
{
  double t1 = sp(p_i, p_a);
  double t2 = sp(p_j, p_a);
  double t3 = sp(p_i, p_j);
  return ((t1 + t2 - t3) / (t1 + t2));
}
// c.f. (5.42) in hep-ph/0201036v1
static inline double Z (const FV& p_i, const FV& p_j, const FV& p_a)
{
  double t1 = sp(p_i, p_a);
  double t2 = sp(p_j, p_a);
  return (t1 / (t1 + t2));
}
// c.f. (5.138) in hep-ph/9605323v3
static inline double x (const FV& p_i, const FV& p_a, const FV& p_b)
{
  double t1 = sp(p_a, p_b);
  double t2 = sp(p_i, p_a);
  double t3 = sp(p_i, p_b);
  return ((t1 - t2 - t3) / t1);
}
// c.f. (5.12) in hep-ph/0201036v1
static inline double y (const FV& p_i, const FV& p_j, const FV& p_k)
{
  double t1 = sp(p_i, p_j);
  double t2 = sp(p_i, p_k);
  double t3 = sp(p_j, p_k);
  return (t1 / (t1 + t2 + t3));
}
// c.f. (5.12) in hep-ph/0201036v1
static inline double z (const FV& p_i, const FV& p_j, const FV& p_k)
{
  double t1 = sp(p_i, p_k);
  double t2 = sp(p_j, p_k);
  return (t1 / (t1 + t2));
}


// = v-tilde_ijk / v_ijk in (5.16)
// c.f. also (5.8), (5.14) in hep-ph/0201036v1
static inline double v (
		 const FV& p_i,
		 const FV& p_j,
		 const FV& p_k,
		 double const& mu2)
{
  double t1 = lambda(0.1e1, mu2, mu2);
  double t2 = sqrt(t1);
  double t3 = y(p_i, p_j, p_k);
  double t4 = 0.1e1 - t3;
  double t10 = pow(0.2e1 * mu2 + (0.1e1 - 0.2e1 * mu2) * t4, 0.2e1);
  double t13 = sqrt(t10 - 0.4e1 * mu2);
  return (t2 * t4 / t13);
}


// diagonal part of the splitting kernel for g->gg from initial state with final spectator
// c.f. (5.85) in hep-ph/0201036v1
// including the propagator 1/(2pa.pi * x_iab) in (5.71)
// but NOT the constant factors -16 pi AlphaS CA * ColorCorr.
static inline double VggIF (FV const& p_a, FV const& p_i, FV const& p_j)
{
  double t1 = sp(p_a, p_i);
  double t3 = X(p_i, p_j, p_a);
  double t6 = Z(p_j, p_i, p_a);
  return (0.1e1 / t1 / t3 * (-0.1e1 + 0.1e1 / (0.2e1 - t3 - t6) + t3 * (0.1e1 - t3)) / 0.2e1);
}
// non-diag. part of the splitting kernel for g->gg from initial state with final spectator
// c.f. (5.85) in hep-ph/0201036v1
// including the propagator 1/(2pa.pi * x_iab) in (5.71)
// but NOT the constant factors -16 pi AlphaS CA * ColorCorr.
static inline double VggIFc (const FV& p_a, const FV& p_i, const FV& p_j)
{
  double t1 = X(p_i, p_j, p_a);
  double t3 = Z(p_i, p_j, p_a);
  double t5 = Z(p_j, p_i, p_a);
  double t7 = sp(p_a, p_i);
  double t9 = t1 * t1;
  double t12 = sp(p_i, p_j);
  return ((0.1e1 - t1) * t3 * t5 / t7 / t9 / t12 / 0.2e1);
}


// diagonal part of the splitting kernel for g->gg from initial state with initial spectator
// c.f. (5.148) in hep-ph/9605323v3
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors -16 pi AlphaS CA * ColorCorr.
static inline double VggII (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = sp(p_a, p_i);
  double t3 = x(p_i, p_a, p_b);
  double t6 = 0.1e1 - t3;
  return (0.1e1 / t1 / t3 * (t3 / t6 + t3 * t6) / 0.2e1);
}
// non-diag. part of the splitting kernel for g->gg from initial state with initial spectator
// c.f. (5.148) in hep-ph/9605323v3
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors -16 pi AlphaS CA * ColorCorr.
static inline double VggIIc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = x(p_i, p_a, p_b);
  double t3 = sp(p_a, p_b);
  double t5 = pow(sp(p_a, p_i),2);
  double t8 = t1 * t1;
  double t13 = sp(p_i, p_b);
  return ((0.1e1 - t1) * t3 / t5 / t8 / t13 / 0.2e1);
}


// diagonal part of the splitting kernel for q->qg from initial state with initial spectator,
// c.f. (5.147) in hep-ph/9605323v3
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors 8 pi AlphaS CF * ColorCorr.
static inline double VqqII (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = sp(p_a, p_i);
  return (-0.1e1 / t1 / 0.2e1);
}
// non-diag. part of the splitting kernel for q->qg from initial state with initial spectator
// c.f. (5.147) in hep-ph/9605323v3
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors 8 pi AlphaS CF * ColorCorr.
static inline double VqqIIc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  return -2.0*VggIIc(p_a,p_i,p_b);
}


// diagonal part of the splitting kernel for q->qg from initial state with final spectator,
// c.f. (5.83) in hep-ph/0201036v1
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors 8 pi AlphaS CF * ColorCorr.
static inline double VqqIF (const FV& p_a, const FV& p_i, const FV& p_b)
{
  return VqqII(p_a,p_i,p_b);
}
// non-diag. part of the splitting kernel for q->qg from initial state with final spectator
// c.f. (5.83) in hep-ph/0201036v1
// including the propagator 1/(2pa.pi * x_iab) in (5.136)
// but NOT the constant factors 8 pi AlphaS CF * ColorCorr.
static inline double VqqIFc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  return -2.0*VggIFc(p_a,p_i,p_b);
}


// diagonal part of the splitting kernel for Q->Qg from final state with final spectator,
// c.f. (5.16) in hep-ph/0201036v1
// including the propagator 1/(2pi.pj) in (5.2)
// but NOT the constant factors -8 pi AlphaS CF * ColorCorr.
static inline double VQgFF (const FV& p_i, const FV& p_j, const FV& p_k)
{
  double t1 = sp(p_i, p_j);
  double t2 = 0.1e1 / t1;
  double t3 = z(p_i, p_j, p_k);
  double t4 = y(p_i, p_j, p_k);
  double t10 = sp(p_i, p_i);
  double t11 = sp(p_i, p_k);
  double t12 = sp(p_j, p_k);
  double t17 = v(p_i, p_j, p_k, t10 / (0.2e1 * t10 + 0.2e1 * t1 + 0.2e1 * t11 + 0.2e1 * t12));
  return (t2 * (0.2e1 / (0.1e1 - t3 * (0.1e1 - t4)) - t17 * (t10 * t2 + t3 + 0.1e1)) / 0.2e1);
}
// diagonal part of the splitting kernel for Q->Qg from final state with initial spectator,
// c.f. (5.50) in hep-ph/0201036v1
// including the propagator 1/(2pi.pj * x_ija) in (5.40)
// but NOT the constant factors -8 pi AlphaS CF * ColorCorr.
static inline double VQgFI (const FV& p_i, const FV& p_j, const FV& p_a)
{
  double t1;
  double t10;
  double t2;
  double t3;
  double t6;
  t1 = sp(p_i, p_j);
  t2 = 0.1e1 / t1;
  t3 = X(p_i, p_j, p_a);
  t6 = Z(p_i, p_j, p_a);
  t10 = sp(p_i, p_i);
  return (t2 / t3 * (-0.1e1 + 0.2e1 / (0.2e1 - t3 - t6) - t6 - t10 * t2) / 0.2e1);
}






// ///////////////////////////////////////////////////////////////////////////////////////////////
// #ifdef WITH_T_SPIN
// ///////////////////////////////////////////////////////////////////////////////////////////////
// // PHI x PHI initial-initial dipole spin correlations
// #include "../amp/real/fullSpin/PHIxPHI_UID_g1E_g2S_t0_tb0.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_UID_g1S_g2E_t0_tb0.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_UID_IM_INTab_g1E_g2S_t0_tb0.cpp"
// #include "../amp/real/fullSpin/PHIxPHI_UID_IM_INTab_g1S_g2E_t0_tb0.cpp"
// // PHI x QCD initial-initial dipole spin correlations
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g2S_t0_tb0.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1S_g2E_t0_tb0.cpp"
// // PHI x QCD initial-final dipole spin correlations
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g20_tS_tb0.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g10_g2E_tS_tb0.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g20_t0_tbS.cpp"
// #include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g10_g2E_t0_tbS.cpp"
// ///////////////////////////////////////////////////////////////////////////////////////////////
// #else
// ///////////////////////////////////////////////////////////////////////////////////////////////
// PHI x PHI initial-initial dipole spin correlations
#include "../amp/real/noSpin/PHIxPHI_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/PHIxPHI_UID_g1S_g2E_t0_tb0.cpp"
#include "../amp/real/noSpin/PHIxPHI_UID_Q_QG_g1E_g2S_t0_tb0.cpp"
// PHI x QCD initial-initial dipole spin correlations
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1S_g2E_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_Q_QG_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_Q_QG_g1E_g20_tS_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_Q_QG_g1E_g20_t0_tbS.cpp"
// PHI x QCD initial-final dipole spin correlations
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g20_tS_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g2E_tS_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g20_t0_tbS.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g2E_t0_tbS.cpp"
// ///////////////////////////////////////////////////////////////////////////////////////////////
// #endif
// ///////////////////////////////////////////////////////////////////////////////////////////////
// //final-final dipoles are proportional to Born
// #include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g20_tE_tbS.cpp"
// #include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g20_tS_tbE.cpp"

// eikonal approximation of the QCD_FSR x (PHI_ISR + PHI_INT) diagrams
#ifdef SGA_WITH_GAUGEVEC
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIISR_INT_SGA_gg.cpp"
#else
#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_PHIISR_INT_SGA_noGaugeVec_gg.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////////////////////
//#endif
///////////////////////////////////////////////////////////////////////////////////////////////

// unintegrated dipoles
double Eval_UID_GG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& norm_factor, // normalization factor for distribution weights
		   DistVec* dist,
		   CutVec*  cuts)
{
  PREFACTORS(hm); // defines ap and hp
  double const& mt2    = hm.mt2();
  double const& mScale = hm.Scale();

  double res = 0.0;

  bool CUTS = (cuts != nullptr);
  bool DIST = (dist != nullptr);

  // 2->3 phase space vectors serve as input for the dipole mapping
  FV const& k1 = ps.k1();
  FV const& k2 = ps.k2();
  FV const& p1 = ps.p1();
  FV const& p2 = ps.p2();
  FV const& p3 = ps.p3();

  // reduced dipole phase space
  // this object is modified in each dipole phase space mapping!!!
  // the initialization is only done once! problems?
  static PS_2_2 ps_red(ps.get_msq(0),ps.get_msq(1),"p1 p2 -> k1 k2 (red.) [static in Eval_UID_GG]");
  static LT lt;
  if (DIST)
    {
      ps_red.P1() = ps.P1();// copy also proton momenta -> needed for lab-frame boost
      ps_red.P2() = ps.P2();// copy also proton momenta -> needed for lab-frame boost
    }
  FV& P1 = ps_red.p1();
  FV& P2 = ps_red.p2();
  FV& K1 = ps_red.k1();
  FV& K2 = ps_red.k2();
  
#ifdef WITH_T_SPIN
  FV const& s1_r = ps.s1_r();
  FV const& s2_r = ps.s2_r();
  FV & S1 = ps_red.s1();
  FV & S2 = ps_red.s2();
#endif
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-INITIAL DIPOLES /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_GG_II(flags))
    {
#endif
      { // DIPOLE-CONFIGURATION: ES00
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================
	double const& bx = x(p3,p1,p2);
	P1 = bx*p1;
	P2 = p2;
	FV Q  = p1 + p2 - p3;
	FV Qt = P1 + p2;
	lt.set_II(Q,Qt);
	// Lorentz transformation of final state vectors
	K1 = k1;
	K2 = k2;
	lt.apply_G(K1);
	lt.apply_G(K2);
#ifdef WITH_T_SPIN
	lt.apply_G(S1);
	lt.apply_G(S2);
#endif
	ps_red.set(bx); // set invariants
	//==============================================================================       		
	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    double PF_VggII  = ap.PREF_UID_CA*VggII(p1,p3,p2);
	    double PF_VggIIc = VggIIc(p1,p3,p2);
#ifdef DEBUG
	    if (flags & F_EVAL_UID_ES00)
	      {
#endif	  
		dip += PF_VggII*Eval_B_2PHIxQCD(ps_red,ap,hp);
		dip += PF_VggIIc*Eval_UID_ES00(ps,ps_red,hm);
#ifdef DEBUG
	      }
	    if (flags & F_EVAL_R_PHIxPHI_ISR)
	      {
#endif	  
		dip += PF_VggII*Eval_B_PHIxPHI_withINT12(ps_red,hm);
		dip += PF_VggIIc*Eval_UID_PHIxPHI_ES00(ps,ps_red,hm);
#ifdef DEBUG
	      }
#endif
	    
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif	  
	    res += dip;
	  }
      }

      { // DIPOLE-CONFIGURATION: SE00
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================
	double const& bx = x(p3,p2,p1);
	P1 = p1;
	P2 = bx*p2;
	FV Q  = p1 + p2 - p3;
	FV Qt = P2 + p1;
	lt.set_II(Q,Qt);
	// Lorentz transformation of final state vectors
	K1 = k1;
	K2 = k2;
	lt.apply_G(K1);
	lt.apply_G(K2);
#ifdef WITH_T_SPIN
        S1 = ps.s1_r();
	S2 = ps.s2_r();
#endif
	ps_red.set(bx);
	//==============================================================================       	

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    double PF_VggII  = ap.PREF_UID_CA*VggII(p2,p3,p1);
	    double PF_VggIIc = VggIIc(p2,p3,p1);
#ifdef DEBUG
	    if (flags & F_EVAL_UID_SE00)
	      {
#endif
		dip += PF_VggII*Eval_B_2PHIxQCD(ps_red,ap,hp);
		dip += PF_VggIIc*Eval_UID_SE00(ps,ps_red,hm);
#ifdef DEBUG
	      }
	    if (flags & F_EVAL_R_PHIxPHI_ISR)
	      {
#endif
		// !!! need to change in case USE2H = true !!!
		dip += PF_VggII*Eval_B_PHIxPHI_withINT12(ps_red,hm);
		dip += PF_VggIIc*Eval_UID_PHIxPHI_SE00(ps,ps_red,hm);
#ifdef DEBUG
	      }
#endif
	
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif	  
	    res += dip;
	  }
      }
#ifdef DEBUG
    }
#endif


  ////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL-FINAL DIPOLES   ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_GG_FF(flags))
    {
#endif
      { // DIPOLE-CONFIGURATION: 00ES
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================
	FV Q   = k1+k2+p3;
	FV kij = k1+p3;
	double mQ2= MSQ(Q);
	P1 = p1;
	P2 = p2;
	// note that this formula is specific for the case that emitter and spectator
	// masses are equal: P_ij^2 = P_k^2 = m_t^2 !!!
	// spectator
	K2 = sqrt(lambda(mQ2,mt2,mt2)/lambda(mQ2,MSQ(kij),mt2))*(k2-(sp(Q,k2)/mQ2)*Q)+0.5*Q;
	// emitter
	K1 = Q-K2;
#ifdef WITH_T_SPIN
	set_spins_in_tt_zmf(K1,K2,S1,S2,s1_r,s2_r);
#endif
	ps_red.set(2.0*(mt2+sp(k1,k2))/ps_red.get_s());
	//==============================================================================

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    double PF_VQgFF = ap.PREF_UID_CF*VQgFF(k1,p3,k2);

#ifdef DEBUG
	    if (flags & F_EVAL_UID_00ES)
	      {
#endif
		dip += PF_VQgFF*Eval_B_2PHIxQCD(ps_red,ap,hp);
#ifdef DEBUG
	      }
	    if (flags & F_EVAL_R_PHIxPHI_FSR)
	      {
#endif
		dip += PF_VQgFF*Eval_B_PHIxPHI_withINT12(ps_red,hm);
#ifdef DEBUG
	      }
#endif
	    
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif	  
	    res += dip;
	  }
      }
  
      { // DIPOLE-CONFIGURATION: 00SE
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================
	FV Q   = k1+k2+p3;
	FV kij = k2+p3;
	double mQ2= MSQ(Q);
	P1 = p1;
	P2 = p2;
	// note that this formula is specific for the case that emitter and spectator
	// masses are equal: P_ij^2 = P_k^2 = m_t^2 !!!
	// spectator
	K1 = sqrt(lambda(mQ2,mt2,mt2)/lambda(mQ2,MSQ(kij),mt2))*(k1-(sp(Q,k1)/mQ2)*Q)+0.5*Q;
	// emitter
	K2 = Q-K1;
#ifdef WITH_T_SPIN
	set_spins_in_tt_zmf(K1,K2,S1,S2,s1_r,s2_r);
#endif
	ps_red.set(2.0*(mt2+sp(k1,k2))/ps_red.get_s());
	//==============================================================================
	
	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);

	    double PF_VQgFF = ap.PREF_UID_CF*VQgFF(k2,p3,k1);

#ifdef DEBUG
	    if (flags & F_EVAL_UID_00SE)
	      {
#endif
		dip += PF_VQgFF*Eval_B_2PHIxQCD(ps_red,ap,hp);
#ifdef DEBUG
	      }
	    if (flags & F_EVAL_R_PHIxPHI_FSR)
	      {
#endif
		dip += PF_VQgFF*Eval_B_PHIxPHI_withINT12(ps_red,hm);
#ifdef DEBUG
	      }
#endif
	    
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
#ifdef DEBUG
    }
#endif

  static const double f = 0.5;
  ////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL-INITIAL DIPOLES ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_GG_FI(flags))
    {
#endif
      // final state emitter
      {
      	// S0E0
      	double dip = 0.0;
      	double const& bx = X(k1,p3,p1);
      	P1 = bx*p1;
      	P2 = p2;
      	K1 = (bx-1.0)*p1 + k1 + p3;
      	K2 = k2;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    //(sp(K1,P2)-sp(K1,P1))/sp(P1,P2);// = -0.5*CA*beta*y
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip += 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k1,p3,p1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// 0SE0
      	double dip = 0.0;
      	double const& bx = X(k1,p3,p2);
      	P2 = bx*p2;
      	P1 = p1;
      	K1 = (bx-1.0)*p2 + k1 + p3;
      	K2 = k2;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip -= 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k1,p3,p2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// S00E
      	double dip = 0.0;
      	double const& bx = X(k2,p3,p1);
      	P1 = bx*p1;
      	P2 = p2;
      	K2 = (bx-1.0)*p1 + k2 + p3;
      	K1 = k1;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip -= 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k2,p3,p1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// 0S0E
      	double dip = 0.0;
      	double const& bx = X(k2,p3,p2);
      	P2 = bx*p2;
      	P1 = p1;
      	K2 = (bx-1.0)*p2 + k2 + p3;
      	K1 = k1;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip += 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k2,p3,p2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
#ifdef DEBUG
    }
#endif

  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-FINAL ///////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_GG_IF(flags))
    {
#endif
      //initial state emitter
      {
      	// E0S0
      	double dip = 0.0;
      	double const& bx = X(p3,k1,p1);
      	P1 = bx*p1;
      	P2 = p2;
      	K1 = (bx-1.0)*p1 + k1 + p3;
      	K2 = k2;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip += 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p1,p3,k1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += VggIFc(p1,p3,k1)*Eval_UID_E0S0(ps,ps_red,hm);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// E00S
      	double dip = 0.0;
      	double const& bx = X(p3,k2,p1);
      	P1 = bx*p1;
      	P2 = p2;
      	K2 = (bx-1.0)*p1 + k2 + p3;
      	K1 = k1;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip -= 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p1,p3,k2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += VggIFc(p1,p3,k2)*Eval_UID_E00S(ps,ps_red,hm);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// 0ES0
      	double dip = 0.0;
      	double const& bx = X(p3,k1,p2);
      	P2 = bx*p2;
      	P1 = p1;
      	K1 = (bx-1.0)*p2 + k1 + p3;
      	K2 = k2;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip -= 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p2,p3,k1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += VggIFc(p2,p3,k1)*Eval_UID_0ES0(ps,ps_red,hm);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
      	// 0E0S
      	double dip = 0.0;
      	double const& bx = X(p3,k2,p2);
      	P2 = bx*p2;
      	P1 = p1;
      	K2 = (bx-1.0)*p2 + k2 + p3;
      	K1 = k1;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {	
	    hm.SetHiggsPrefactors(ps_red.get_s(),1);
	    dip += 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p2,p3,k2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += VggIFc(p2,p3,k2)*Eval_UID_0E0S(ps,ps_red,hm);
	    dip *= f;

	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
#ifdef DEBUG
    }
#endif

  
#ifndef WITH_T_SPIN
#ifdef DEBUG
  if (EVAL_SGA_GG(flags))
    {
#endif
      //the real FSRxISR + FSRxINT contributions are subtracted by their SGA
      if (EvalCuts(cuts,&ps,nullptr))
	{
	  hm.SetHiggsPrefactors(2.0*(mt2+sp(k1,k2)),1);
#ifdef SGA_WITH_GAUGEVEC
	  static FV gv;
	  gv[0] = +p3[0];
	  gv[1] = -p3[1];
	  gv[2] = -p3[2];
	  gv[3] = -p3[3];
	  double sga = Eval_R_FSR_ISR_SGA(ps,ap,hp,gv);
#else
	  double sga = Eval_R_FSR_ISR_SGA(ps,ap,hp);
#endif
	  // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	  if (DIST)
	    {
	      // in case of the SGA collect distributions from 2->3 phase space
	      ps.FillDistributions(*dist,H_NLO_PHI_R,sga*norm_factor,mScale);
	    }
	  ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(sga);
#endif
	  res+= sga;
	}
#ifdef DEBUG
    }
#endif // DEBUG
#endif // WITH_T_SPIN

  return res;
}








double Eval_UID_QG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& norm_factor, // normalization factor for distribution weights
		   DistVec* dist,
		   CutVec*  cuts )
{
  PREFACTORS(hm); // defines ap and hp
  using namespace Constants;
  double const& mScale = hm.Scale();
  
  double res = 0.0;

  bool DIST = (dist != nullptr);

  // 2->3 phase space vectors serve as input for the dipole mapping
  FV const& k1 = ps.k1();
  FV const& k2 = ps.k2();
  FV const& p1 = ps.p1();
  FV const& p2 = ps.p2();
  FV const& p3 = ps.p3();

  // reduced dipole phase space
  // this object is modified in each dipole phase space mapping!!!
  // the initialization is only done once! problems?
  static PS_2_2 ps_red(ps.get_msq(0),ps.get_msq(1),"p1 p2 -> k1 k2 (red.) [static in Eval_UID_QG]");
  static LT lt;
  if (DIST)
    {
      ps_red.P1() = ps.P1();// copy also proton momenta -> needed for lab-frame boost
      ps_red.P2() = ps.P2();// copy also proton momenta -> needed for lab-frame boost
    }
  FV& P1 = ps_red.p1();
  FV& P2 = ps_red.p2();
  FV& K1 = ps_red.k1();
  FV& K2 = ps_red.k2();
  
#ifdef WITH_T_SPIN
  FV const& s1 = ps.s1();
  FV const& s2 = ps.s2();
  FV & S1 = ps_red.s1();
  FV & S2 = ps_red.s2();
#endif
  ////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-INITIAL DIPOLES /////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_QG_II(flags))
    { // DIPOLE-CONFIGURATION: ES00
#endif
      double dip = 0.0;
      //==============================================================================
      // dipole phase space mapping
      //==============================================================================
      double bx = x(p3,p1,p2);
      P1 = bx*p1;
      P2 = p2;
      FV Q  = p1 + p2 - p3;
      FV Qt = P1 + p2;
      lt.set_II(Q,Qt);
      // Lorentz transformation of final state vectors
      lt.apply_G_cpy(k1,K1);
      lt.apply_G_cpy(k2,K2);
#ifdef WITH_T_SPIN
      lt.apply_G_cpy(s1,S1);
      lt.apply_G_cpy(s2,S2);
#endif
      ps_red.set(bx);


      if (EvalCuts(cuts,&ps_red,nullptr))
	{
	  hm.SetHiggsPrefactors(ps_red.get_s(),1);
	  FV PI = p3-sp(p1,p3)/sp(p1,p2)*p2;
    	
	  double Vqq_diag = 8.0*Pi*hm.AlphaS()*CF/CA*VqqII(p1,p3,p2);
	  double Vqq_nond = 8.0*Pi*hm.AlphaS()*CF/CA*VqqIIc(p1,p3,p2);

	  if (flags & F_EVAL_R_PHIxQCD_QG)
	    {
	      dip += Eval_UID_Q_QG_ES00(ps,ps_red,hm,PI,Vqq_diag,Vqq_nond);
	    }
	  if (flags & F_EVAL_R_PHIxPHI_QG)
	    {
	      dip += Eval_UID_PHIxPHI_Q_QG_ES00(ps,ps_red,hm,PI,Vqq_diag,Vqq_nond);
	    }

	  // works
	  // double Vgg_diag = -16.0*Pi*hm.AlphaS()*VggII(p1,p3,p2);
	  // double Vgg_nond = -16.0*Pi*hm.AlphaS()*VggIIc(p1,p3,p2);

	  // if (flags & F_EVAL_R_ISR_ISR)
	  // 	{
	  // 	  dip += Eval_UID_Q_QG_ES00(ps,ps_red,hm,PI,Vgg_diag,Vgg_nond);
	  // 	}
	  // if (flags & F_EVAL_R_PHIxPHI_ISR)
	  // 	{
	  // 	  dip += Eval_UID_PHIxPHI_Q_QG_ES00(ps,ps_red,hm,PI,Vgg_diag,Vgg_nond);
	  // 	}      

	  // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	  if (DIST)
	    {
	      ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	    }
	  ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	  CHECKNAN(dip);
#endif	  
	  res += dip;
	}
#ifdef DEBUG
    }
#endif
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-FINAL DIPOLES ///////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
  if (EVAL_UID_QG_IF(flags))
    {
#endif
      {
	// DIPOLE-CONFIGURATION: E0S0
	double dip = 0.0;
	double bx  = X(p3,k1,p1);
	P1 = bx*p1;
	P2 = p2;
	K1 = (bx-1.0)*p1 + k1 + p3;
	K2 = k2;
	ps_red.set(bx);
    
	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),true);	
	    double Vqq_diag = 8.0*Pi*hm.AlphaS()*CF/CA*VqqIF(p1,p3,k1);
	    double Vqq_nond = 8.0*Pi*hm.AlphaS()*CF/CA*VqqIFc(p1,p3,k1);

	    FV PI = 1.0/Z(p3,k1,p1)*p3-1.0/Z(k1,p3,p1)*k1;
	
	    if (flags & F_EVAL_R_PHIxQCD_QG)
	      {
		dip += Eval_UID_Q_QG_E0S0(ps,ps_red,hm,PI,Vqq_diag,Vqq_nond);
	      }
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
      {
	// DIPOLE-CONFIGURATION: E00S
	double dip = 0.0;
	double bx  = X(p3,k2,p1);
	P1 = bx*p1;
	P2 = p2;
	K2 = (bx-1.0)*p1 + k2 + p3;
	K1 = k1;
	ps_red.set(bx);

	if (EvalCuts(cuts,&ps_red,nullptr))
	  {
	    hm.SetHiggsPrefactors(ps_red.get_s(),true);	
	    double Vqq_diag = 8.0*Pi*hm.AlphaS()*CF/CA*VqqIF(p1,p3,k2);
	    double Vqq_nond = 8.0*Pi*hm.AlphaS()*CF/CA*VqqIFc(p1,p3,k2);

	    FV PI = 1.0/Z(p3,k2,p1)*p3-1.0/Z(k2,p3,p1)*k2;
	
	    if (flags & F_EVAL_R_PHIxQCD_QG)
	      {
		dip += Eval_UID_Q_QG_E00S(ps,ps_red,hm,PI,Vqq_diag,Vqq_nond);
	      }
	    // DISTRIBUTIONS ///////////////////////////////////////////////////////////
	    if (DIST)
	      {
		ps_red.FillDistributions(*dist,H_NLO_PHI_R,dip*norm_factor,mScale);
	      }
	    ////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	    CHECKNAN(dip);
#endif
	    res += dip;
	  }
      }
#ifdef DEBUG
    }
#endif
  return res;
}













