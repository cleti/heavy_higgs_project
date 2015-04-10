
#include "../../../../inc/Functions_Shared.h"


double Eval_B_2PHIxQCD
(
 FV const& p1,
 FV const& p2,
 FV const& k1,
 FV const& k2
 )
{
  const double& s = 2.0*sp(p1,p2);
  const double& t11 = 2.0*sp(p1,k1); 
  const double& t12 = 2.0*sp(p1,k2); 

  c_double t2 = DenS(s,mH,GammaH);
  double t1;
  double t14;
  double t23;
  double t3;
  double t4;
  t1 = FH0*RE(t2);
  t3 = t12 / 0.2e1 + t11 / 0.2e1;
  t4 = t3 * t3;
  t14 = 0.1e1 / t12 / t11;
  t23 = FA0*RE(t2);
  return(-0.128e3 * t1 * t4 * (-0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t14 * PREF_B_At + 0.256e3 * t23 * t4 * t3 * t14 * PREF_B_Bt);
}

