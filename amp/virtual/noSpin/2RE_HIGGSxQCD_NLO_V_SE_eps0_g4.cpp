
#include "AMP_HEADER.h"

#include <math.h>

double Eval_V_SE (AMP_ARGS)
{

  AMP_DEFINITIONS
  HP_REFS_PHIxQCD(hp)
  AP_REFS_B(ap)
  AP_REFS_V(ap)

  double t1;
  double t11;
  double t13;
  double t15;
  double t16;
  double t19;
  double t2;
  double t21;
  double t22;
  double t24;
  double t25;
  double t27;
  double t28;
  double t29;
  double t3;
  double t30;
  double t31;
  double t32;
  double t33;
  double t34;
  double t35;
  double t36;
  double t39;
  double t4;
  double t40;
  double t42;
  double t43;
  double t44;
  double t48;
  double t50;
  double t53;
  double t54;
  double t56;
  double t6;
  double t60;
  double t65;
  double t68;
  double t7;
  double t73;
  double t74;
  double t75;
  double t8;
  double t80;
  double t81;
  double t82;
  double t85;
  double t9;
  double t93;
  double t94;
  double t95;
  double t96;
  double t97;
  t1 = beta * beta;
  t2 = t1 * t1;
  t3 = s * s;
  t4 = t3 * s;
  t6 = y * y;
  t7 = t6 * t6;
  t8 = t2 * t4 * t7;
  t9 = t1 * beta;
  t11 = t6 * y;
  t13 = 0.3e1 * t9 * t4 * t11;
  t15 = t9 * t3 * t11;
  t16 = 0.8e1 * t15;
  t19 = 0.5e1 * t1 * t4 * t6;
  t21 = t1 * t3 * t6;
  t22 = 0.8e1 * t21;
  t24 = t1 * s * t6;
  t25 = 0.16e2 * t24;
  t27 = 0.5e1 * beta_y * t4;
  t28 = beta_y * t3;
  t29 = 0.20e2 * t28;
  t30 = 0.2e1 * t4;
  t31 = beta_y * s;
  t32 = 0.32e2 * t31;
  t33 = 0.20e2 * t3;
  t34 = 0.64e2 * s;
  t35 = t8 - t13 - t16 + t19 + t22 - t25 - t27 + t29 + t30 - t32 - t33 + t34 - 0.64e2;
  t36 = 0.1e1 / s;
  t39 = pow(beta_y - 0.1e1, 0.2e1);
  t40 = 0.1e1 / t39;
  t42 = 0.1e1 / (t31 - s + 0.2e1);
  t43 = t40 * t42;
  t44 = RE(I2_T11_0_MT2_MU2_0);
  t48 = t8 + t13 + t16 + t19 + t22 - t25 + t27 - t29 + t30 + t32 - t33 + t34 - 0.64e2;
  t50 = 0.1e1 / (t31 + s - 0.2e1);
  t53 = pow(beta_y + 0.1e1, 0.2e1);
  t54 = 0.1e1 / t53;
  t56 = RE(I2_T12_0_MT2_MU2_0);
  t60 = t24 - s + 0.4e1;
  t65 = t6 * t1;
  t68 = 0.8e1 * s;
  t73 = t40 * t36;
  t74 = RE(I1_MT2_MU2_0);
  t75 = t50 * t74;
  t80 = t2 * s * t7;
  t81 = 0.2e1 * t24;
  t82 = 0.4e1 * t65;
  t85 = 0.3e1 * (RE(I2_MT2_0_MT2_MU2_0)-2.0) + 0.5e1;
  // RE(I2_MT2_0_MT2_MU2_0)-2.0 = log(MUR^2/mt^2)
  t93 = 0.3e1 * t15;
  t94 = 0.7e1 * t21;
  t95 = 0.4e1 * t24;
  t96 = 0.5e1 * t28;
  t97 = 0.4e1 * t31;
  return(At_fH_re * (0.8e1 * t35 * t36 * t43 * t44 - 0.8e1 * t48 * t50 * t36 * t54 * t56 - 0.32e2 * t60 * (t2 * t3 * t7 - 0.3e1 * t21 + 0.12e2 * t24 - 0.2e1 * t3 - 0.8e1 * t65 + t68 - 0.8e1) * t54 * t42 * t73 * t75 + 0.64e2 * (t80 - t81 + t82 - s + 0.4e1) * t85 * t73 * t54) * PREF_V_CF + Bt_fA_re * (-0.16e2 * (t93 - t94 + t95 + t96 + t97 - t3 - t68 + 0.16e2) * t40 * t42 * t44 - 0.16e2 * (t93 + t94 - t95 + t96 + t97 + t3 + t68 - 0.16e2) * t50 * t54 * t56 + 0.64e2 * (-0.3e1 * t24 + 0.2e1 * t65 - s + 0.2e1) * t60 * t54 * t43 * t75 + 0.32e2 * (t80 - t81 + t82 + s + 0.4e1) * t85 * t40 * t54) * PREF_V_CF);
}
