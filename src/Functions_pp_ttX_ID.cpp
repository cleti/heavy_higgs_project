
#include "../inc/Functions_pp_ttX_ID.h"



#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\






// regularized Altarelli-Parisi splitting function g->gg, c.f. (5.89) in hep-ph/0201036v1
inline double P_gg_reg(double const& x)
{
  double t1 = 0.1e1 - x;
  return Constants::CA*2.0*( t1 / x - 1.0 +  t1 * x); // changed 29.4.
}

// regularized  Altarelli-Parisi splitting function q->qg, c.f. (5.89) in hep-ph/0201036v1
inline double P_qg_reg(double const& x) 
{
  return Constants::CF*(1.0+std::pow(1.0-x,2))/x;
}
// eps^1 coefficient in the d-dimensional q->qg splitting function, c.f. (5.93) in hep-ph/0201036v1
inline double P1_qg(double const& x) 
{
  return Constants::CF*(x);
}

#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_GG_DELTA_I_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_GG_DELTA_PK_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_GG_CONT_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_GG_DIST_E_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_QG_CONT_P_eps0_g4.cpp"



// evaluate the finite part of the integrated dipoles: x-independent part of I,P,K operators
double Eval_ID_GG(
	       const PS_2_2& ps,
	       HiggsModel& hm,
	       const ulong& flags)
{
  if (flags & (F_EVAL_D_GG_DELTA | F_EVAL_D_GG_END))
    {
      PREFACTORS(hm); // defines ap and hp
      double res = 0.0;

      //using namespace RunParameters;
      const double& MUR2   = hm.MUR2();
      const double& MUF2   = hm.MUF2();
      const double& AlphaS = hm.AlphaS();
      const double& s12_1  = ps.get_s();


      // these have to be adjusted before calls to Eval_??? whenever S changes!!!
      hm.SetHiggsPrefactors(s12_1,1);
      double B_phi = 0.0;
      if (flags & F_EVAL_B_PHIxPHI) 
	{
	  B_phi += Eval_B_PHIxPHI(ps,ap,hp);
	}
      double B_int = 0.0;
      if (flags & F_EVAL_B_PHIxQCD) 
	{
	  B_int += Eval_B_2PHIxQCD(ps,ap,hp);
	}
  
      // colour connected Born terms: only interference terms PHIxQCD are non-zero
      double const& beta_1 = ps.get_beta();
      double const& y_1    = ps.get_y();
      double BC_int        = -0.5*beta_1*y_1*B_int;

      // x-independent part (terms proportional to delta(1-x))
      double const& S12_1  = 2.0*sp(ps.k1(),ps.k2());
      double const& t11_1  = ps.get_t11();
      double const& t12_1  = ps.get_t12();


      if (flags & F_EVAL_D_GG_DELTA) 
	{
	  res += Eval_DIP_GG_DELTA_I (s12_1,
	  			   S12_1,
	  			   beta_1,
	  			   t11_1,
	  			   t12_1,
	  			   B_phi+B_int,
	  			   BC_int,
	  			   MUR2,
	  			   AlphaS);
	  // res += Eval_DIP_GG_DELTA_PK(s12_1,
	  // 			   t11_1,
	  // 			   t12_1,
	  // 			   B_phi+B_int,
	  // 			   2.0*B_phi+B_int,
	  // 			   BC_int,
	  // 			   MUF2,MUR2,
	  // 			   AlphaS);
	}
      // distribution end-point part (check tech. cut also here)
      double const& x      = ps.get_x();
      if (flags & F_EVAL_D_GG_END && x<Cuts::IDIP_X_CUT) 
	{
	  res += Eval_DIP_GG_DIST_E(s12_1,
				    t11_1,
				    t12_1,
				    B_phi+B_int,    // II and FF colour correlated Born MEs
				    // these are proportional to the uncorrelated Born MEs (factors -CA/CF)
				    2.0*B_phi+B_int,// this is the Born contribution that is multiplied with the K-bar terms
				    // the IF colour correlated Born contribution for PHIxPHI is zero, but the K-bar terms
				    // of the IF dipoles do not involve colour operators thus PHIxPHI contributes!
				    // further the factor 1/2 has to be compensated for the PHIxPHI part
				    // this factor was introduced because half of the IF phiXqcd dipoles were neglected
				    BC_int, // IF and FI colour correlated Born MEs
				    // only the PHIxQCD interference terms contribute here
				    x,
				    MUF2,
				    AlphaS);
	}
      return res;
    }
  else
    {
      return 0.0;
    }
}








// evaluate the finite part of the integrated dipoles: continuous, x-dependent contribution of P and K operators
double Eval_ID_X_GG(
		 const PS_2_2& ps_x,
		 HiggsModel& hm,
		 const ulong& flags)
{

  if ( flags & (F_EVAL_D_GG_CONT) )
    {
      double const& x    = ps_x.get_x();
      // check technical cut
      if (x>Cuts::IDIP_X_CUT) return 0.0;

      PREFACTORS(hm); // defines ap and hp
      
      double const& MUF2   = hm.MUF2();
      double const& AlphaS = hm.AlphaS();

      // c.m.e. on boosted phase space
      double const& s12_x = ps_x.get_s();
      // original partonic c.m.e.
      double s12_1        = s12_x / x;

      // invariants t_aj = 2 p_a.k_j, where p_a is the original initial state 
      // momentum and k_j comes from the boosted phase space ps_x
      // c.f. text below (6.38) in hep-ph/0201036v1
      double t11_x = ps_x.get_t11() / x; 
      double t12_x = ps_x.get_t12() / x; 

      // these have to be adjusted before calls to Eval_??? whenever S changes!!!
      hm.SetHiggsPrefactors(s12_x,1);
      double B_phi_x = 0.0;
      if (flags & F_EVAL_B_PHIxPHI) 
	{
	  B_phi_x += Eval_B_PHIxPHI(ps_x,ap,hp);
	}
#ifdef DEBUG
      CHECKNA(B_phi_x);
#endif
      
      double B_int_x = 0.0;
      double BC_int_x= 0.0;
      if (flags & F_EVAL_B_PHIxQCD) 
	{
	  B_int_x += Eval_B_2PHIxQCD(ps_x,ap,hp);
	  // colour correlated Born MEs, only PHIxQCD contributes
	  BC_int_x+= -0.5*ps_x.get_beta_y()*B_int_x;
	}
#ifdef DEBUG
      CHECKNA(B_int_x);
#endif

      return Eval_DIP_GG_CONT(s12_1,
			      t11_x,
			      t12_x,
			      B_phi_x+B_int_x,
			      2.0*B_phi_x+B_int_x,
			      BC_int_x,
			      x,
			      MUF2,
			      AlphaS);
    }
  else
    {
      return 0.0;
    }
}


double Eval_ID_X_QG(
		    const PS_2_2& ps_x,
		    HiggsModel& hm,
		    const ulong& flags)
{
  if ( flags & (F_EVAL_D_QG_CONT) )
    {
      double const& x    = ps_x.get_x();
      // check technical cut
      if (x>Cuts::IDIP_X_CUT) return 0.0;

      PREFACTORS(hm); // defines ap and hp
      
      double const& MUF2   = hm.MUF2();
      double const& AlphaS = hm.AlphaS();

      // c.m.e. on boosted phase space
      double const& s12_x = ps_x.get_s();
      // original partonic c.m.e.
      double s12_1        = s12_x / x;

      // invariants t_aj = 2 p_a.k_j, where p_a is the original initial state 
      // momentum and k_j comes from the boosted phase space ps_x
      // c.f. text below (6.38) in hep-ph/0201036v1
      double t11_x = ps_x.get_t11() / x; 
      double t12_x = ps_x.get_t12() / x; 

      // these have to be adjusted before calls to Eval_??? whenever S changes!!!
      hm.SetHiggsPrefactors(s12_x,1);

      // I HAVE CHANGED THIS RECENTLY, CHECK IF SOMETHING IS WRONG
      // PHI^2 + PHIxQCD Born interference terms
      double B_phi_x = 0.0;
      if (flags & F_EVAL_B_PHIxPHI) 
	{
	  B_phi_x += Eval_B_PHIxPHI(ps_x,ap,hp);
	}
#ifdef DEBUG
      CHECKNA(B_phi_x);
#endif
      
      double B_int_x = 0.0;
      double BC_int_x= 0.0;
      if (flags & F_EVAL_B_PHIxQCD) 
	{
	  B_int_x += Eval_B_2PHIxQCD(ps_x,ap,hp);
	  // colour correlated Born MEs, only PHIxQCD contributes
	  BC_int_x+= -0.5*ps_x.get_beta_y()*B_int_x;
	}
#ifdef DEBUG
      CHECKNA(B_int_x);
#endif     

      // the dipoles associated with the q->qg splitting only have a continuum part
      return Eval_DIP_QG_CONT(s12_1,
			      t11_x,
			      t12_x,
			      B_phi_x+B_int_x,
			      BC_int_x,
			      x,
			      MUF2,
			      AlphaS);
    }
  else
    {
      return 0.0;
    }
}
