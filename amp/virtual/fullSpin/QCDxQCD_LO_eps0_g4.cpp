
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t10;
  double t102;
  double t12;
  double t13;
  double t15;
  double t16;
  double t18;
  double t19;
  double t20;
  double t22;
  double t23;
  double t29;
  double t31;
  double t35;
  double t39;
  double t4;
  double t41;
  double t42;
  double t44;
  double t47;
  double t48;
  double t49;
  double t50;
  double t57;
  double t58;
  double t6;
  double t61;
  double t62;
  double t68;
  double t69;
  double t71;
  double t72;
  double t8;
  double t84;
  double t87;
  double t9;
  double t96;
  double t99;
  t1 = beta * beta;
  t4 = y * y;
  t6 = t1 * beta * s * t4 * y;
  t8 = t1 * s * t4;
  t9 = beta_y * s;
  t10 = 0.4e1 * beta_y;
  t12 = s * s;
  t13 = 0.1e1 / t12;
  t15 = beta_y - 0.1e1;
  t16 = 0.1e1 / t15;
  t18 = beta_y + 0.1e1;
  t19 = 0.1e1 / t18;
  t20 = sp(s1, p1);
  t22 = sp(k1, s2);
  t23 = t19 * t20 * t22;
  t29 = sp(s2, p1);
  t31 = sp(s1, k2);
  t35 = 0.1e1 / s;
  t39 = t1 * t1;
  t41 = t4 * t4;
  t42 = t39 * s * t41;
  t44 = t1 * t4;
  t47 = t35 * (t42 - 0.2e1 * t8 + 0.8e1 * t44 + s);
  t48 = t16 * t19;
  t49 = sp(s1, s2);
  t50 = t48 * t49;
  t57 = 0.4e1 * t44;
  t58 = 0.2e1 * beta_y;
  t61 = t15 * t15;
  t62 = 0.1e1 / t61;
  t68 = t18 * t18;
  t69 = 0.1e1 / t68;
  t71 = t16 * t29;
  t72 = t71 * t31;
  t84 = t8 - s + 0.4e1;
  t87 = t62 * t69;
  t96 = t39 * t12 * t41;
  t99 = t1 * t12 * t4;
  t102 = 0.4e1 * s;
  return((-0.64e2 * (t6 + t8 - t9 + t10 - s) * t13 * t16 * t23 + 0.64e2 * (t6 - t8 - t9 + t10 + s) * t13 * t16 * t19 * t29 * t31 + 0.128e3 * t35 * t29 * t20 - 0.16e2 * t47 * t50 + 0.16e2 * t47 * t48) * PREF_B_QCDxQCD_CA + (0.256e3 * (t42 - t8 + t57 - t58 + 0.2e1) * t13 * t62 * t23 - 0.256e3 * (t42 - t8 + t57 + t58 + 0.2e1) * t13 * t69 * t72 - 0.512e3 * (t42 - t8 + 0.2e1 * t44 + 0.2e1) * t13 * t62 * t69 * t29 * t20 + 0.64e2 * (t42 - t8 + t57 + 0.4e1) * t84 * t13 * t87 * t49 - 0.64e2 * (t1 * t12 * t39 * t4 * t41 - t102 - t12 + 0.8e1 * t42 + 0.16e2 * t44 - 0.4e1 * t8 - t96 + t99 + 0.16e2) * t13 * t87) * PREF_B_QCDxQCD_CF + (-0.256e3 * (t6 + t8 - t9 + t10 - s + 0.2e1) * t13 * t19 * t16 * t20 * t22 + 0.256e3 * (t6 - t8 - t9 + t10 + s - 0.2e1) * t13 * t19 * t72 + 0.512e3 * (t8 - s + 0.2e1) * t13 * t19 * t71 * t20 - 0.64e2 * (t96 - 0.2e1 * t99 + 0.8e1 * t8 + t12 - t102 + 0.16e2) * t13 * t50 + 0.64e2 * t84 * (t8 + 0.4e1) * t13 * t19 * t16) * PREF_B_QCDxQCD_CFCA);
}
