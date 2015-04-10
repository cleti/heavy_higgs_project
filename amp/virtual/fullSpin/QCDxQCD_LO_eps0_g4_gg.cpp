
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_GG (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t10;
  double t11;
  double t12;
  double t14;
  double t15;
  double t2;
  double t20;
  double t21;
  double t24;
  double t26;
  double t27;
  double t3;
  double t32;
  double t33;
  double t34;
  double t35;
  double t37;
  double t38;
  double t39;
  double t4;
  double t40;
  double t42;
  double t46;
  double t48;
  double t5;
  double t57;
  double t59;
  double t6;
  double t7;
  double t73;
  t1 = beta * beta;
  t2 = y - 0.1e1;
  t3 = t1 * t2;
  t4 = y + 0.1e1;
  t5 = beta_y - 0.1e1;
  t6 = t5 * t5;
  t7 = 0.1e1 / t6;
  t10 = beta_y + 0.1e1;
  t11 = 0.1e1 / t10;
  t12 = 0.1e1 / s;
  t14 = sp(s2, p2);
  t15 = sp(s1, p1);
  t20 = t10 * t10;
  t21 = 0.1e1 / t20;
  t24 = 0.1e1 / t5;
  t26 = sp(s2, p1);
  t27 = sp(s1, p2);
  t32 = t1 * t1;
  t33 = y * y;
  t34 = t33 * t33;
  t35 = t32 * t34;
  t37 = 0.2e1 * t32 * t33;
  t38 = 0.2e1 * t32;
  t39 = 0.2e1 * t1;
  t40 = t35 - t37 + t38 - t39 + 0.1e1;
  t42 = sp(s1, s2);
  t46 = t1 * t33;
  t48 = t35 - t37 + t38 + 0.2e1 * t46 - t39 - 0.1e1;
  t57 = t46 + 0.1e1;
  t59 = t2 * t4 * t1 * t57;
  t73 = t7 * t21;
  return((-0.128e3 * t11 * t12 * t14 * t15 * t3 * t4 * t7 + 0.128e3 * t12 * t21 * t24 * t26 * t27 * t3 * t4 - 0.32e2 * t21 * t40 * t42 * t7 + 0.32e2 * t21 * t48 * t7) * PREF_B_QCDxQCD / CA + (0.64e2 * t11 * t12 * t14 * t15 * t59 * t7 - 0.64e2 * t12 * t21 * t24 * t26 * t27 * t59 + 0.16e2 * t40 * t42 * t57 * t73 - 0.16e2 * t48 * t57 * t73) * PREF_B_QCDxQCD_CA);
}
