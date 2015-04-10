
#include "../inc/Functions_pp_HX.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define AMP_PP_HX_HEADER			\
  using namespace Constants;			\
  using namespace AmpPrefactors;		\
  using namespace RunParameters;		\
  using namespace HiggsBosons;			\
  double const FH02 = FH_eff_1*FH_eff_1;        \
  double const FA02 = FA_eff_1*FA_eff_1;	\


double Eval_B(
	      FV const& p1,
	      FV const& p2
	      )
{
  AMP_PP_HX_HEADER

  double t5 = PREF_B_PHIxPHI/CA;
  return (t5 * (std::norm(F_ggH1_s) + std::norm(F_ggH1_p)));
}



double Eval_B_eff(
		  FV const& p1,
		  FV const& p2
		  )
{
  AMP_PP_HX_HEADER

  double s = 2.0*sp(p1,p2);
  double t10 = s * s;
  double t6 = 16.0 * PREF_B_PHIxPHI/CA;
  return (t6 * (0.4e1 * FA02 + FH02) * t10);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Eval_V(
	      FV const& p1,
	      FV const& p2
	      )
{
  AMP_PP_HX_HEADER


  double s = 2.0*sp(p1,p2);  
  double t1 = s * s;
  // double t4 = log(MUR2 / s);
  double t5 = 0.0;//t4 * t4; // cancels exactly against term in I-op.
  double t6 = CA * t5;
  double t10 = CA * Pi2;
  return 0.8e1 * t1 * (0.16e2 * CA * FA02 + 0.4e1 * t10 * FA02 - 0.4e1 * t6 * FA02 + t10 * FH02 - t6 * FH02 + 0.11e2 * FH02) / CA * PREF_V_PHI;
}


static double P_gg_reg(double const& x)
{
  static double x_t   = 0.0;
  static double res_t = 0.0;
  if (x!=x_t)
    {
      x_t = x;
      double t1 = 0.1e1 - x;
      res_t = Constants::CA*(0.2e1 / x * t1 - 0.2e1 + 0.2e1 * t1 * x);
    }
  return res_t;
}





#include "../amp/ggHX/2RE_HIGGSxQCD_NLO_DIP_DELTA_I_eps0_g4.cpp"
#include "../amp/ggHX/2RE_HIGGSxQCD_NLO_DIP_DELTA_PK_eps0_g4.cpp"
#include "../amp/ggHX/2RE_HIGGSxQCD_NLO_DIP_CONT_eps0_g4.cpp"
#include "../amp/ggHX/2RE_HIGGSxQCD_NLO_DIP_DIST_E_eps0_g4.cpp"

// evaluate the finite part of the integrated dipoles: x-independent part of I,P,K operators
double Eval_ID(
	       FV const& p1,
	       FV const& p2,
	       double const& x
	       )
{
  using namespace RunParameters;
  double res = 0.0;

  double B_phi = Eval_B_eff(p1,p2);
  double const& s12_1  = 2.0*sp(p1,p2);
  res += Eval_DIP_DELTA_I (s12_1, // invairant in initial-initial dipoles
  			   0.0,   // invariant in final-final dipoles
  			   0.0,   // beta factor for massive final-final dipoles
  			   0.0,   // invariant t11 (final-initial)
  			   0.0,   // invariant t12 (final-initial)
  			   B_phi, // uncorrelated Born 
  			   0.0,   // color correlated Born
  			   MUR2); // ren./fact. scales
  res += Eval_DIP_DELTA_PK(s12_1, // invariant s12 (initial-initial) 
  			   0.0,   // invariant t12 (initial-final)	 
  			   0.0,   // invariant t12 (initial-final)   
  			   B_phi, // uncorrelated Born	       
  			   0.0,   // uncorrelated Born in IF dipoles 
  			   0.0,   // colour correlated Born          
  			   MUF2,MUR2);// ren./fact. scales

  // distribution end-point part
  if (x<Cuts::IDIP_X_CUT) 
    {
      res -= Eval_DIP_DIST_E(s12_1, // invariant s12 (initial-initial) 
      			     0.0,   // invariant t12 (initial-final)   	 
      			     0.0,   // invariant t12 (initial-final)   
      			     B_phi, // uncorrelated Born	       
      			     0.0,   // uncorrelated Born in IF dipoles 
      			     0.0,   // colour correlated Born          
      			     x,     // boost factor		       
      			     MUF2); // ren./fact. scales       
      
    }
#ifdef DEBUG
  CHECKNA(res);
#endif
  return res;
}



// evaluate the finite part of the integrated dipoles: continuous, x-dependent contribution of P and K operators
double Eval_ID_X(
		 FV const& p1,
		 FV const& p2,
		 double const& x
		 )
{
  using namespace RunParameters;
  double res = 0.0;
  // this contributes only if the rescaled c.m.e. is still above threshold (boost_initial_state() checks this)
  if ( x<Cuts::IDIP_X_CUT )
    {
      double rx = sqrt(x);
      double B_phi_x = Eval_B_eff(rx*p1,rx*p2);
      double s12_1   = 2.0*sp(p1,p2);
      res += Eval_DIP_CONT(s12_1, // invariant s12 (initial-initial) 
			   0.0,	    // invariant t12 (initial-final)   
			   0.0,	    // invariant t12 (initial-final)   
			   B_phi_x, // uncorrelated Born	       
			   0.0,	    // uncorrelated Born in IF dipoles 
			   0.0,	    // colour correlated Born          
			   x,	    // boost factor		       
			   MUF2);   // ren./fact. scales               
    }
#ifdef DEBUG
  CHECKNA(res);
#endif
  return res;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Eval_R_gg (
	       FV const& p1,
	       FV const& p2,
	       FV const& p3
	       )
{
  AMP_PP_HX_HEADER

  double t4 = sp(p1, p2);
  double t5 = sp(p1, p3);
  double t6 = sp(p2, p3);
  double t9 = 0.16e2 * pow(t4 - t5 - t6, 0.4e1);
  double t10 = t4 * t4;
  double t11 = t10 * t10;
  double t13 = t5 * t5;
  double t14 = t13 * t13;
  double t16 = t6 * t6;
  double t17 = t16 * t16;
  return 0.16e2 * PREF_R_PHI * (4.0 * FA02 + FH02) * (t9 + 0.16e2 * t11 + 0.16e2 * t14 + 0.16e2 * t17) / (t4 * t5 * t6);
}






static double x (const FV& p_i, const FV& p_a, const FV& p_b)
{
  double t1 = sp(p_a, p_b);
  double t2 = sp(p_i, p_a);
  double t3 = sp(p_i, p_b);
  return ((t1 - t2 - t3) / t1);
}

static double VggII (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = sp(p_a, p_i);
  double t3 = x(p_i, p_a, p_b);
  double t6 = 0.1e1 - t3;
  return (0.1e1 / t1 / t3 * (t3 / t6 + t3 * t6) / 0.2e1);
}

static double VggIIc (const FV& p_a, const FV& p_i, const FV& p_b)
{
  double t1 = x(p_i, p_a, p_b);
  double t3 = sp(p_a, p_b);
  double t5 = pow(sp(p_a, p_i),2);
  double t8 = t1 * t1;
  double t13 = sp(p_i, p_b);
  return ((0.1e1 - t1) * t3 / t5 / t8 / t13 / 0.2e1);
}




double Eval_UID_II_C_ES (
			FV const& p1,
			FV const& p2,
			FV const& p3,
			FV const& P1,
			FV const& P2
			)
{
  AMP_PP_HX_HEADER

  double t1 = sp(P1,p2);
  double t2 = sp(P1,p3);
  double t3 = sp(p2,p3);
  return 1024.0 * PREF_R_PHI * (FH02 + 4.0 * FA02) * t1 * t2 * t3;
}





double Eval_UID (
		 FV const& p1,
		 FV const& p2,
		 FV const& p3
		 )
{
  double res = 0.0;
  static FV P1;
  static FV P2;

  // p1: emitter, p2: spectator
  {
    P1 = x(p3,p1,p2)*p1;
    P2 = p2;
    res += AmpPrefactors::PREF_UID_CA*VggII(p1,p3,p2)*Eval_B_eff(P1,P2);
    res += Eval_UID_II_C_ES(p1,p2,p3,P1,p2)*VggIIc(p1,p3,p2);
  }
  // p2: emitter, p1: spectator
  {
    P1 = p1;
    P2 = x(p3,p2,p1)*p2;
    res += AmpPrefactors::PREF_UID_CA*VggII(p2,p3,p1)*Eval_B_eff(P1,P2);
    res += Eval_UID_II_C_ES(p2,p1,p3,P2,p1)*VggIIc(p2,p3,p1);
  }
#ifdef DEBUG
  CHECKNA(res);
#endif  
  return res;
}




// probably wrong
double Eval_R_qq(
		 FV const& p1,
		 FV const& p2,
		 FV const& p3
		 )
{
  AMP_PP_HX_HEADER
    
  double t16;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t5 = PREF_V_CT;
  t6 = sp(p3, p2);
  t7 = t6 * t6;
  t8 = sp(p1, p3);
  t9 = t8 * t8;
  t16 = sp(p1, p2);
  return (0.128e3 * t5 * (t7 + t9) * (0.4e1 * FA02 + FH02) / t16);
}



double Eval_sigma_R_qq_fin(
		     double const& z
		     )
{
  AMP_PP_HX_HEADER
   
  double t1 = pow(1-z,3);
  double t18 = (0.4e1 * FA02 + FH02);
  return (AlphaS2*AlphaS/Pi2*144.0/486.0 * t18  * t1);
}

// phase space, spin/color average and flux factors already included!!!
// note:  z = mH2/s 
double Eval_sigma_R_qg_fin(
		     double const& z
		     )
{
  AMP_PP_HX_HEADER
  
  double P_qg = CF*(1.0+std::pow(1.0-z,2))/z;
  double t18 = (0.4e1 * FA02 + FH02);
  return (AlphaS2*AlphaS/(Pi2*4.0)) * t18 * ((z-1.0)*(7.0-3.0*z)/3.0 + 0.5*z*P_qg*(1.0+std::log(M2_1/MUF2*std::pow(1.0-z,2)/z)) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


