
#include "../inc/Functions_pp_ttX_UID.h"




#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\



static double X (const FV&p_i, const FV&p_j, const FV&p_a)
{
  double t1 = sp(p_i, p_a);
  double t2 = sp(p_j, p_a);
  double t3 = sp(p_i, p_j);
  return ((t1 + t2 - t3) / (t1 + t2));
}

static double Z (const FV& p_i, const FV& p_j, const FV& p_a)
{
  double t1 = sp(p_i, p_a);
  double t2 = sp(p_j, p_a);
  return (t1 / (t1 + t2));
}

static double x (const FV& p_i, const FV& p_a, const FV& p_b)
{
  double t1 = sp(p_a, p_b);
  double t2 = sp(p_i, p_a);
  double t3 = sp(p_i, p_b);
  return ((t1 - t2 - t3) / t1);
}

static double y (const FV& p_i, const FV& p_j, const FV& p_k)
{
  double t1 = sp(p_i, p_j);
  double t2 = sp(p_i, p_k);
  double t3 = sp(p_j, p_k);
  return (t1 / (t1 + t2 + t3));
}

static double z (const FV& p_i, const FV& p_j, const FV& p_k)
{
  double t1 = sp(p_i, p_k);
  double t2 = sp(p_j, p_k);
  return (t1 / (t1 + t2));
}



static double v (
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


// splitting kernel for g->gg from initial state with final spectator
static double VggIF (FV const& p_a, FV const& p_i, FV const& p_j)
{
  double t1 = sp(p_a, p_i);
  double t3 = X(p_i, p_j, p_a);
  double t6 = Z(p_j, p_i, p_a);
  return (0.1e1 / t1 / t3 * (-0.1e1 + 0.1e1 / (0.2e1 - t3 - t6) + t3 * (0.1e1 - t3)) / 0.2e1);
}
// splitting kernel for g->gg from initial state with final spectator (spin-correlations)
static double VggIFc (const FV& p_a, const FV& p_i, const FV& p_j)
{
  double t1 = X(p_i, p_j, p_a);
  double t3 = Z(p_i, p_j, p_a);
  double t5 = Z(p_j, p_i, p_a);
  double t7 = sp(p_a, p_i);
  double t9 = t1 * t1;
  double t12 = sp(p_i, p_j);
  return ((0.1e1 - t1) * t3 * t5 / t7 / t9 / t12 / 0.2e1);
}
// splitting kernel for g->gg from initial state with initial spectator
static double VggII (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = sp(p_a, p_i);
  double t3 = x(p_i, p_a, p_b);
  double t6 = 0.1e1 - t3;
  return (0.1e1 / t1 / t3 * (t3 / t6 + t3 * t6) / 0.2e1);
}
// splitting kernel for g->gg from initial state with initial spectator (spin-correlations)
static double VggIIc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = x(p_i, p_a, p_b);
  double t3 = sp(p_a, p_b);
  double t5 = pow(sp(p_a, p_i),2);
  double t8 = t1 * t1;
  double t13 = sp(p_i, p_b);
  return ((0.1e1 - t1) * t3 / t5 / t8 / t13 / 0.2e1);
}


// splitting kernel for q->qg from initial state with initial spectator
static double VqqII (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = sp(p_a, p_i);
  return(0.1e1 / t1 / 0.2e1);
}
// splitting kernel for q->qg from initial state with initial spectator (spin-correlations)
static double VqqIIc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  // double t1;
  // double t10;
  // double t13;
  // double t3;
  // double t5;
  // double t8;
  // t1 = x(p_i, p_a, p_b);
  // t3 = sp(p_a, p_b);
  // t5 = sp(p_a, p_i);
  // t8 = t1 * t1;
  // t10 = sp(p_i, p_a);
  // t13 = sp(p_i, p_b);
  // return((0.1e1 - t1) * t3 / t5 / t8 / t10 / t13);
  // only valid in 4 spacetime dimensions !!!
  return 2.0*VggIIc(p_a,p_i,p_b);
}

static double VQgFF (const FV& p_i, const FV& p_j, const FV& p_k)
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

static double VQgFI (const FV& p_i, const FV& p_j, const FV& p_a)
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






///////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
///////////////////////////////////////////////////////////////////////////////////////////////
// PHI x PHI initial-initial dipole spin correlations
#include "../amp/real/fullSpin/PHIxPHI_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/fullSpin/PHIxPHI_UID_g1S_g2E_t0_tb0.cpp"
#include "../amp/real/fullSpin/PHIxPHI_UID_IM_INTab_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/fullSpin/PHIxPHI_UID_IM_INTab_g1S_g2E_t0_tb0.cpp"
// PHI x QCD initial-initial dipole spin correlations
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1S_g2E_t0_tb0.cpp"
// PHI x QCD initial-final dipole spin correlations
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g20_tS_tb0.cpp"
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g10_g2E_tS_tb0.cpp"
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g1E_g20_t0_tbS.cpp"
#include "../amp/real/fullSpin/2RE_PHIxQCD_UID_g10_g2E_t0_tbS.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////
#else
///////////////////////////////////////////////////////////////////////////////////////////////
// PHI x PHI initial-initial dipole spin correlations
#include "../amp/real/noSpin/PHIxPHI_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/PHIxPHI_UID_g1S_g2E_t0_tb0.cpp"
#include "../amp/real/noSpin/PHIxPHI_UID_Q_QG_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/PHIxPHI_UID_Q_QG_g1S_g2E_t0_tb0.cpp"
// PHI x QCD initial-initial dipole spin correlations
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1S_g2E_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_Q_QG_g1E_g2S_t0_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_Q_QG_g1S_g2E_t0_tb0.cpp"
// PHI x QCD initial-final dipole spin correlations
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g20_tS_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g2E_tS_tb0.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g1E_g20_t0_tbS.cpp"
#include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g2E_t0_tbS.cpp"
///////////////////////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////////////////////
// //final-final dipoles are proportional to Born
// #include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g20_tE_tbS.cpp"
// #include "../amp/real/noSpin/2RE_PHIxQCD_UID_g10_g20_tS_tbE.cpp"

// eikonal approximation of the QCD_FSR x (PHI_ISR + PHI_INT) diagrams
//#include "../amp/real/noSpin/2RE_PHIxQCD_NLO_R_QCDFSR_HISR_gg_soft.cpp"
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
		   std::vector<HistArray*>* dist )
{
  PREFACTORS(hm); // defines ap and hp
  double const& mt2    = hm.mt2();
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
  static PS_2_2 ps_red(ps.get_msq(0),ps.get_msq(1),"p1 p2 -> k1 k2 (red.) [static in Eval_UID_GG]");
  static LT lt;
  ps_red.P1() = ps.P1();// copy also proton momenta -> needed for lab-frame boost
  ps_red.P2() = ps.P2();// copy also proton momenta -> needed for lab-frame boost
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
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-INITIAL DIPOLES //////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  if (EVAL_UID_GG_II(flags))
    {
      { // DIPOLE-CONFIGURATION: ES00
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================       
	P1 = x(p3,p1,p2)*p1;
	P2 = p2;
	FV Q  = p1 + p2 - p3;
	FV Qt = P1 + p2;
	lt.set_II(Q,Qt);
	// Lorentz transformation of final state vectors
	K1 = k1;
	K2 = k2;
	lt.apply_G(K1);
	lt.apply_G(K2);
	////////////////////////////
	// boost to parton c.m.f. //
	////////////////////////////
	// lt.set_boost(P1+P2);
	// lt.apply(P1);
	// lt.apply(P2);
	// lt.apply(K1);
	// lt.apply(K2);
	////////////////////////////
	ps_red.set();
#ifdef WITH_T_SPIN
	S1 = ps.s1_r();
	lt.set_boost(K1,1);
	lt.apply_G(S1);
	S2 = ps.s2_r();
	lt.set_boost(K2,1);
	lt.apply_G(S2);
#endif
	////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT NOTE:
	// the 4-vectors of the reduced phase space are not defined in the parton c.m.f.
	// instead they are boosted along the z-axis
	// for distributions of Lorentz non-invariants this has to be taken into account
	//==============================================================================	
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
	double PF_VggII  = ap.PREF_UID_CA*VggII(p1,p3,p2);
	double PF_VggIIc = VggIIc(p1,p3,p2);
	if (flags & F_EVAL_UID_ES00)
	  {
	    dip += PF_VggII*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += PF_VggIIc*Eval_UID_ES00(ps,ps_red,ap,hp);
	  }
	if (flags & F_EVAL_R_PHIxPHI_ISR)
	  {
	    dip += PF_VggII*Eval_B_PHIxPHI_withINT12(ps_red,hm);
	    dip += PF_VggIIc*Eval_UID_PHIxPHI_ES00(ps,ps_red,ap,hp);
#ifdef WITH_T_SPIN
	    if (TwoHDM>0)
	      {
		dip += PF_VggIIc*Eval_UID_PHIxPHI_IM_INTab_SE00(ps,ps_red,ap,hp);
	      }
#endif
	  }
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-initial 2->2 impulskonfig. (ES00): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
	
	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // PRINT_4VEC(K1);
	    // PRINT_4VEC(K2);
	    // PRINT_4VEC(P1);
	    // PRINT_4VEC(P2);
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    // lt.apply(P1);
	    // lt.apply(P2);
	    // PRINT_4VEC(K1);
	    // PRINT_4VEC(K2);
	    // PRINT_4VEC(P1);
	    // PRINT_4VEC(P2);
	    // exit(1);
	    ////////////////////////////
	    // ps_red.P1() = 1.0/x(p3,p1,p2)*ps.P1();
	    // ps_red.P2() = ps.P2();
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif	  
	res += dip;
      }

      { // DIPOLE-CONFIGURATION: SE00
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================       
	P1 = p1;
	P2 = x(p3,p2,p1)*p2;
	FV Q  = p1 + p2 - p3;
	FV Qt = P2 + p1;
	lt.set_II(Q,Qt);
	// Lorentz transformation of final state vectors
	K1 = k1;
	K2 = k2;
	lt.apply_G(K1);
	lt.apply_G(K2);
	////////////////////////////
	// boost to parton c.m.f. //
	////////////////////////////
	// lt.set_boost(P1+P2);
	// lt.apply(P1);
	// lt.apply(P2);
	// lt.apply(K1);
	// lt.apply(K2);
	////////////////////////////
	ps_red.set();
#ifdef WITH_T_SPIN
        S1 = ps.s1_r();
	// lt.set_boost(K1,1);
	// lt.apply_G(S1);
	S2 = ps.s2_r();
	// lt.set_boost(K2,1);
	// lt.apply_G(S2);
#endif
	//==============================================================================       	
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
	double PF_VggII  = ap.PREF_UID_CA*VggII(p2,p3,p1);
	double PF_VggIIc = VggIIc(p2,p3,p1);
	if (flags & F_EVAL_UID_SE00)
	  {
	    dip += PF_VggII*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += PF_VggIIc*Eval_UID_SE00(ps,ps_red,ap,hp);
	  }
	if (flags & F_EVAL_R_PHIxPHI_ISR)
	  {
	    // !!! need to change in case USE2H = true !!!
	    dip += PF_VggII*Eval_B_PHIxPHI_withINT12(ps_red,hm);
	    dip += PF_VggIIc*Eval_UID_PHIxPHI_SE00(ps,ps_red,ap,hp);
#ifdef WITH_T_SPIN
	    if (TwoHDM>0)
	      {
		dip += PF_VggIIc*Eval_UID_PHIxPHI_IM_INTab_SE00(ps,ps_red,ap,hp);
	      }
#endif
	  }
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-initial 2->2 impulskonfig. (SE00): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
	
	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    // ps_red.P1() = ps.P1();
	    // ps_red.P2() = 1.0/x(p3,p2,p1)*ps.P2();
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif	  
	res += dip;
      }
    }


  /////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL-FINAL DIPOLES   ////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  if (EVAL_UID_GG_FF(flags))
    {
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
	// spectator
	// note that this formula is specific for the case that emitter and spectator
	// masses are equal: P_ij^2 = P_k^2 = m_t^2 !!!
	K1 = sqrt(lambda(mQ2,mt2,mt2)/lambda(mQ2,MSQ(kij),mt2))*(k2-(sp(Q,k2)/mQ2)*Q)+0.5*Q;
	// emitter
	K2 = Q-K1;	
	// strange: according to Catani/Seymour Eq. (5.9) the mapping should be the other way round (k1<->k2)
	ps_red.set();
#ifdef WITH_T_SPIN
	set_spins_in_tt_zmf(K1,K2,S1,S2,s1_r,s2_r);
#endif
	//==============================================================================

	double PF_VQgFF = ap.PREF_UID_CF*VQgFF(k1,p3,k2);
	hm.SetHiggsPrefactors(ps.get_s(),1);
	if (flags & F_EVAL_UID_00ES)
	  {
	    dip += PF_VQgFF*Eval_B_2PHIxQCD(ps_red,ap,hp);
	  }
	if (flags & F_EVAL_R_PHIxPHI_FSR)
	  {
	    // !!! need to change in case USE2H = true !!!
	    dip += PF_VQgFF*Eval_B_PHIxPHI_withINT12(ps_red,hm);
	  }
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-final 2->2 impulskonfig. (00ES): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
	
	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    // ps_red.print();
	    // ps.print();
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif	  
	res += dip;
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
	// spectator
	// note that this formula is specific for the case that emitter and spectator
	// masses are equal: P_ij^2 = P_k^2 = m_t^2 !!!
	K2 = sqrt(lambda(mQ2,mt2,mt2)/lambda(mQ2,MSQ(kij),mt2))*(k1-sp(Q,k1)/mQ2*Q)+0.5*Q;
	// emitter
	K1 = Q-K2;
	ps_red.set();
#ifdef WITH_T_SPIN
	set_spins_in_tt_zmf(K1,K2,S1,S2,s1_r,s2_r);
#endif	
	//==============================================================================

	double PF_VQgFF = ap.PREF_UID_CF*VQgFF(k2,p3,k1);
	hm.SetHiggsPrefactors(ps.get_s(),1);
	if (flags & F_EVAL_UID_00SE)
	  {
	    dip += PF_VQgFF*Eval_B_2PHIxQCD(ps_red,ap,hp);
	  }
	if (flags & F_EVAL_R_PHIxPHI_FSR)
	  {
	    // !!! need to change in case USE2H = true !!!
	    dip += PF_VQgFF*Eval_B_PHIxPHI_withINT12(ps_red,hm);
	  }

#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-final 2->2 impulskonfig. (00SE): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif

	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
	////////////////////////////////////////////////////////////////////////////
	res += dip;
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      }
    }

  /////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-FINAL and FINAL-INITIAL DIPOLES //////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  if (flags & (F_EVAL_R_ISR_FSR))
    {

      static double f = 0.5;
      // final state emitter
      {
      	// S0E0
      	double dip = 0.0;
      	double const& x = X(k1,p3,p1);
      	P1 = x*p1;
      	P2 = p2;
      	K1 = (x-1.0)*p1 + k1 + p3;
      	K2 = k2;
	ps_red.set();
	//(sp(K1,P2)-sp(K1,P1))/sp(P1,P2);// = -0.5*CA*beta*y
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip += 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k1,p3,p1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip *= f;
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-initial 2->2 impulskonfig. (S0E0): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// 0SE0
      	double dip = 0.0;
      	double const& x = X(k1,p3,p2);
      	P2 = x*p2;
      	P1 = p1;
      	K1 = (x-1.0)*p2 + k1 + p3;
      	K2 = k2;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip -= 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k1,p3,p2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip *= f;

#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-initial 2->2 impulskonfig. (0SE0): ";
	std::cout << std::endl << "----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// S00E
      	double dip = 0.0;
      	double const& x = X(k2,p3,p1);
      	P1 = x*p1;
      	P2 = p2;
      	K2 = (x-1.0)*p1 + k2 + p3;
      	K1 = k1;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip -= 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k2,p3,p1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip *= f;

#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-initial 2->2 impulskonfig. (S00E): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// 0S0E
      	double dip = 0.0;
      	double const& x = X(k2,p3,p2);
      	P2 = x*p2;
      	P1 = p1;
      	K2 = (x-1.0)*p2 + k2 + p3;
      	K1 = k1;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip += 0.25*ps_red.get_beta_y()*ap.PREF_UID_CA*VQgFI(k2,p3,p2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip *= f;
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " final-initial 2->2 impulskonfig. (0S0E): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }

      //initial state emitter
      {
      	// E0S0
      	double dip = 0.0;
      	double const& x = X(p3,k1,p1);
      	P1 = x*p1;
      	P2 = p2;
      	K1 = (x-1.0)*p1 + k1 + p3;
      	K2 = k2;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip += 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p1,p3,k1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip += VggIFc(p1,p3,k1)*Eval_UID_E0S0(ps,ps_red,ap,hp);
      	dip *= f;

#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-final 2->2 impulskonfig. (E0S0): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// E00S
      	double dip = 0.0;
      	double const& x = X(p3,k2,p1);
      	P1 = x*p1;
      	P2 = p2;
      	K2 = (x-1.0)*p1 + k2 + p3;
      	K1 = k1;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip -= 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p1,p3,k2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip += VggIFc(p1,p3,k2)*Eval_UID_E00S(ps,ps_red,ap,hp);
      	dip *= f;
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-final 2->2 impulskonfig. (E00S): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// 0ES0
      	double dip = 0.0;
      	double const& x = X(p3,k1,p2);
      	P2 = x*p2;
      	P1 = p1;
      	K1 = (x-1.0)*p2 + k1 + p3;
      	K2 = k2;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip -= 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p2,p3,k1)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip += VggIFc(p2,p3,k1)*Eval_UID_0ES0(ps,ps_red,ap,hp);
      	dip *= f;
		
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-final 2->2 impulskonfig. (0ES0): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }
      {
      	// 0E0S
      	double dip = 0.0;
      	double const& x = X(p3,k2,p2);
      	P2 = x*p2;
      	P1 = p1;
      	K2 = (x-1.0)*p2 + k2 + p3;
      	K1 = k1;
	ps_red.set();
	hm.SetHiggsPrefactors(ps_red.get_s(),1);
      	dip += 0.5*ps_red.get_beta_y()*ap.PREF_UID_CA*VggIF(p2,p3,k2)*Eval_B_2PHIxQCD(ps_red,ap,hp);
      	dip += VggIFc(p2,p3,k2)*Eval_UID_0E0S(ps,ps_red,ap,hp);
      	dip *= f;
			
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-final 2->2 impulskonfig. (0E0S): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
      	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif
      	res += dip;
      }    
    }
#ifndef WITH_T_SPIN
  if (flags & (F_EVAL_R_FSR_ISR))
    {
      //the real FSRxISR + FSRxINT contributions are subtracted by their SGA
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
	    ps.FillDistributions(*dist,H_NLO_R,sga*norm_factor,mScale);
	  }
      	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(sga);
#endif
      	res+= sga;
      }
    }
#endif
  return res;
}








double Eval_UID_QG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& norm_factor, // normalization factor for distribution weights
		   std::vector<HistArray*>* dist )
{
  PREFACTORS(hm); // defines ap and hp
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
  ps_red.P1() = ps.P1();// copy also proton momenta -> needed for lab-frame boost
  ps_red.P2() = ps.P2();// copy also proton momenta -> needed for lab-frame boost
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

  /////////////////////////////////////////////////////////////////////////////////////////////
  // INITIAL-INITIAL DIPOLES //////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  if (EVAL_UID_QG_II(flags))
    {
      { // DIPOLE-CONFIGURATION: ES00
	double dip = 0.0;
	//==============================================================================
	// dipole phase space mapping
	//==============================================================================       
	P1 = x(p3,p1,p2)*p1;
	P2 = p2;
	FV Q  = p1 + p2 - p3;
	FV Qt = P1 + p2;
	lt.set_II(Q,Qt);
	// Lorentz transformation of final state vectors
	K1 = k1;
	K2 = k2;
	lt.apply_G(K1);
	lt.apply_G(K2);
	////////////////////////////
	// boost to parton c.m.f. //
	////////////////////////////
	// lt.set_boost(P1+P2);
	// lt.apply(P1);
	// lt.apply(P2);
	// lt.apply(K1);
	// lt.apply(K2);
	////////////////////////////
	ps_red.set();
#ifdef WITH_T_SPIN
	S1 = ps.s1_r();
	lt.set_boost(K1,1);
	lt.apply_G(S1);
	S2 = ps.s2_r();
	lt.set_boost(K2,1);
	lt.apply_G(S2);
#endif
	////////////////////////////////////////////////////////////////////////////////
	// IMPORTANT NOTE:
	// the 4-vectors of the reduced phase space are not defined in the parton c.m.f.
	// instead they are boosted along the z-axis
	// for distributions of Lorentz non-invariants this has to be taken into account
	//==============================================================================	

	hm.SetHiggsPrefactors(ps_red.get_s(),1);	
	double PF_VqqII  = ap.PREF_UID_TF*VqqII(p1,p3,p2);
	double PF_VqqIIc = VqqIIc(p1,p3,p2);

	if (flags & F_EVAL_R_PHIxQCD_QG)
	  {
	    dip += PF_VqqII*Eval_B_2PHIxQCD(ps_red,ap,hp);
	    dip += PF_VqqIIc*Eval_UID_Q_QG_ES00(ps,ps_red,ap,hp);
	  }
	if (flags & F_EVAL_R_PHIxPHI_QG)
	  {
	    dip += PF_VqqII*Eval_B_PHIxPHI_withINT12(ps_red,hm);
	    dip += PF_VqqIIc*Eval_UID_PHIxPHI_Q_QG_ES00(ps,ps_red,ap,hp);
#ifdef WITH_T_SPIN
	    // if (TwoHDM>0)
	    //   {
	    // 	dip += PF_VggIIc*Eval_UID_PHIxPHI_IM_INTab_SE00(ps,ps_red);
	    //   }
#endif	    
	  }
	
#ifdef DUMP_DIPOLE_PS
	std::cout << std::endl << "-----------------------------------------------------------------";
	std::cout << std::endl << " initial-initial 2->2 impulskonfig. (ES00): ";
	std::cout << std::endl << "-----------------------------------------------------------------";
	PRINT_4VEC(mScale*P1);
	PRINT_4VEC(mScale*P2);
	PRINT_4VEC(mScale*K1);
	PRINT_4VEC(mScale*K2);
	std::cout << "-----------------------------------------------------------------" << std::endl;
#endif
	
	// DISTRIBUTIONS ///////////////////////////////////////////////////////////
	if (DIST)
	  {
	    ////////////////////////////
	    // boost to parton c.m.f. //
	    ////////////////////////////
	    // lt.set_boost(P1+P2);
	    // lt.apply(K1);
	    // lt.apply(K2);
	    ////////////////////////////
	    ps_red.FillDistributions(*dist,H_NLO_R,dip*norm_factor,mScale);
	  }
	////////////////////////////////////////////////////////////////////////////
#ifdef DEBUG
	CHECKNAN(dip);
#endif	  
	res += dip;
      }
    }
  return res;
}












