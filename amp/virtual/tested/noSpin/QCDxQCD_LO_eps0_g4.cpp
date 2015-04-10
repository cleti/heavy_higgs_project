
#include "../../../../inc/Functions_Shared.h"

double Eval_B_QCDxQCD (
		       PS_2_2 const& ps
		       )
{

#include "EVAL_V_PS_REFS"

  double t1;
  double t14;
  double t17;
  double t2;
  double t3;
  t1 = beta_y * beta_y;
  t2 = t1 * t1;
  t3 = s * s;
  t14 = pow(t1 - 0.1e1, 2);
  t17 = t14 * t3;
  return(-0.16e2 * (0.8e1 * s * t1 + t2 * t3 - 0.8e1 * s - t3 + 0.32e2) * ( ( t1 - 0.1e1 ) * PREF_B_QCDxQCD_CA + 0.4e1 * PREF_B_QCDxQCD_CF) / t17 );
  
  // double t1;
  // double t10;
  // double t11;
  // double t13;
  // double t19;
  // double t21;
  // double t3;
  // double t4;
  // t1 = s * s;
  // t3 = beta_y * beta_y;
  // t4 = t3 * t3;
  // t10 = 0.1e1 / t1 * (0.8e1 * s * t3 + t1 * t4 - 0.8e1 * s - t1 + 0.32e2);
  // t11 = beta_y + 0.1e1;
  // t13 = beta_y - 0.1e1;
  // t19 = t13 * t13;
  // t21 = t11 * t11;
  // return(-0.16e2 * t10 / t11 / t13 * PREF_B_QCDxQCD_CA - 0.64e2 * t10 / t19 / t21 * PREF_B_QCDxQCD_CF);
}
