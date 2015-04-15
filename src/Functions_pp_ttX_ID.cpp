
#include "../inc/Functions_pp_ttX_ID.h"



#define PREFACTORS(HM)						\
  const AmplitudePrefactors& ap = HM.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = HM.GetHiggsPrefactors();	\





static double P_gg_reg(double const& x)
{
  double t1 = 0.1e1 - x;
  return Constants::CA*(0.2e1 / x * t1 - 0.2e1 + 0.2e1 * t1 * x);
}


inline double P_qg_reg(double const& x) // = P_qg
{
  double t1 = 0.1e1 - x;
  return Constants::CF*(1+t1*t1)/x;
}


#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_DELTA_I_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_DELTA_PK_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_CONT_eps0_g4.cpp"
#include "../amp/virtual/dipole/2RE_HIGGSxQCD_NLO_DIP_DIST_E_eps0_g4.cpp"


inline double Eval_DIP_CONT_QG(
			double const& s12_1,
			double const& B,
			double const& x,
			double const& MUF2,
			double const& AS)
{
  using namespace Constants;
  return -AS/(TwoPi*2.0*CF)*P_qg_reg(x)*std::log(MUF2/(s12_1*x*(1.0-x)))*B;
}




// evaluate the finite part of the integrated dipoles: x-independent part of I,P,K operators
double Eval_ID(
	       const PS_2_2& ps,
	       HiggsModel& hm,
	       const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp
  
  // FV const& p1 = ps.p1();
  // FV const& p2 = ps.p2();
  FV const& k1 = ps.k1();
  FV const& k2 = ps.k2();
  //using namespace RunParameters;
  const double& MUR2 = hm.MUR2();
  const double& MUF2 = hm.MUF2();
  const double& AS   = hm.AlphaS();
  const double& s12_1 = ps.get_s();
  
  double res = 0.0;

  // these have to be adjusted before calls to Eval_??? whenever S changes!!!
  //HiggsBosons::ResetHiggsPrefactors(s12_1);
  hm.SetHiggsPrefactors(s12_1,1);
  // I HAVE CHANGED THIS RECENTLY, CHECK IF SOMETHING IS WRONG
  // PHIxQCD Born interference terms
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
  double const& S12_1  = 2.0*sp(k1,k2);
  double const& t11_1  = ps.get_t11();
  double const& t12_1  = ps.get_t12();
  if (flags & F_EVAL_D_DELTA) 
    {
      res += Eval_DIP_DELTA_I (s12_1,
      			       S12_1,
      			       beta_1,
      			       t11_1,
      			       t12_1,
      			       B_phi+B_int,
      			       BC_int,
      			       MUR2,
			       AS);
      res += Eval_DIP_DELTA_PK(s12_1,
      			       t11_1,
      			       t12_1,
      			       B_phi+B_int,
      			       2.0*B_phi+B_int,// this is the Born contribution that is multiplied with the K-bar terms
      			       BC_int,
      			       MUF2,MUR2,
			       AS);
    }
  // distribution end-point part
  double const& x      = ps.get_x();
  if (flags & F_EVAL_D_END && x<Cuts::IDIP_X_CUT) 
    {
      res += Eval_DIP_DIST_E(s12_1,
			     t11_1,
			     t12_1,
			     B_phi+B_int,
			     2.0*B_phi+B_int,// this is the Born contribution that is multiplied with the K-bar terms
			     BC_int,
			     x,
			     MUF2,
			     AS);
    }

  // the IF dipoles give a P+K contribution for PHI^2, although the colour correlations are zero
  // there is a part in the K operator that is proportional to the uncorrelated Born MEs
  // there, the factor 1/2 has to becompensated for the PHI^2 part
  // this factor was introduced because half of the dipoles contructed from the interference terms only are neglected
  
  return res;
}








// evaluate the finite part of the integrated dipoles: continuous, x-dependent contribution of P and K operators
double Eval_ID_X(
		 const PS_2_2& ps_x,
		 HiggsModel& hm,
		 const ulong& flags)
{
  PREFACTORS(hm); // defines ap and hp
  
  double const& x    = ps_x.get_x();
  //using namespace RunParameters;
  double const& MUF2 = hm.MUF2();
  const double& AS   = hm.AlphaS();
  
  if (flags & (F_EVAL_D_CONT | F_EVAL_D_QG) && x<Cuts::IDIP_X_CUT)
    {
      double ret = 0.0;
      // FV const& p1   = ps_1.p1();
      // FV const& p2   = ps_1.p2();
      // FV const& k1   = ps_1.k1();
      // FV const& k2   = ps_1.k2();

      double const& s_x = ps_x.get_s();
      double s12_1      = s_x / x;
      
      FV const& p1_x = ps_x.p1();
      FV const& p2_x = ps_x.p2();
      FV const& k1_x = ps_x.k1();
      FV const& k2_x = ps_x.k2();

      // invariants of ps_x 4 vectors
      // these are defined in a reference frame that is boosted in z-direction with respect to the z.m.f. of ps !
      double s12_x = 2.0*sp(p1_x,p2_x) / x; // = 2.0*sp(p1,p2_b)  equivalent to the 
      double t11_x = 2.0*sp(p1_x,k1_x) / x; // = 2.0*sp(p1,k1_b)  definition of the invariants
      double t12_x = 2.0*sp(p1_x,k2_x) / x; // = 2.0*sp(p1,k2_b)  in the Catani/Seymour paper


      // these have to be adjusted before calls to Eval_??? whenever S changes!!!
      //HiggsBosons::ResetHiggsPrefactors(s_x);
      hm.SetHiggsPrefactors(s_x,1);
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
      if (flags & F_EVAL_B_PHIxQCD) 
	{
	  B_int_x += Eval_B_2PHIxQCD(ps_x,ap,hp);
	}
#ifdef DEBUG
      CHECKNA(B_int_x);
#endif     

      // colour connected Born terms: only contribution from interference terms PHIxQCD are non-zero
      double BC_int_x= -0.5*(t12_x-t11_x)/s12_x*B_int_x;

      ret += Eval_DIP_CONT(s12_1,
			   t11_x,
			   t12_x,
			   B_phi_x+B_int_x,
			   2.0*B_phi_x+B_int_x,// this is the Born contribution that is multiplied with the K-bar terms
			   BC_int_x,
			   x,
			   MUF2,
			   AS);

      if (flags & F_EVAL_D_QG) ret += Eval_DIP_CONT_QG(s12_1,B_phi_x+B_int_x,x,MUF2,AS);

      return ret;
    }
  else
    {
      return 0.0;
    }
}


