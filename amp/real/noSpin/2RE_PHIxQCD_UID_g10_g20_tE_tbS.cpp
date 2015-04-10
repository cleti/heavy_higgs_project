
#include "UID_HEADER.h"

double Eval_UID_00ES (UID_ARGS)
{

UID_DEFINITIONS


  double t10;
  double t11;
  double t12;
  double t14;
  double t15;
  double t16;
  double t20;
  double t22;
  double t3;
  double t33;
  double t34;
  double t4;
  double t44;
  double t45;
  double t5;
  double t7;
  double t8;
  t3 = VQgFF(k1, p3, k2);
  t4 = sp(K2, p2);
  t5 = t3 * t4;
  t7 = sp(p1, p2);
  t8 = t3 * t7;
  t10 = sp(K1, p1);
  t11 = 0.1e1 / t10;
  t12 = t11 * t3;
  t14 = sp(K1, p2);
  t15 = 0.1e1 / t14;
  t16 = t15 * t3;
  t20 = sp(K1, K2);
  t22 = t3 * t14;
  t33 = sp(K2, p1);
  t34 = t3 * t33;
  t44 = EPS_(K1, K2, p1, p2);
  t45 = t3 * t44;
  return(At_fH_re * (0.256e3 * t5 + 0.512e3 * t8 + (-0.256e3 * t12 * t7 - 0.256e3 * t16 * t7) * t20 - 0.512e3 * t22 + (-0.256e3 * t5 * t7 + 0.256e3 * t8) * t15 + (0.256e3 * t16 * t4 - 0.512e3 * t3) * t10 + 0.256e3 * t34 + (0.256e3 * t22 * t33 - 0.256e3 * t34 * t7 + 0.256e3 * t8) * t11) * PREF_R_CF + At_fA_re * (-0.512e3 * t11 * t45 + 0.512e3 * t15 * t45) * PREF_R_CF + Bt_fA_re * (0.512e3 * t12 * t33 * t7 + 0.512e3 * t16 * t4 * t7 + 0.1024e4 * t8) * PREF_R_CF);
}
