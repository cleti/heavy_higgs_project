
#include "AMP_HEADER.h"

double Eval_V_4G (AMP_ARGS)
{

AMP_DEFINITIONS


  double t10;
  double t103;
  double t110;
  double t111;
  double t114;
  double t116;
  double t12;
  double t128;
  double t13;
  double t15;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t22;
  double t23;
  double t26;
  double t27;
  double t28;
  double t3;
  double t30;
  double t32;
  double t34;
  double t38;
  double t39;
  double t4;
  double t41;
  double t43;
  double t51;
  double t55;
  double t60;
  double t61;
  double t63;
  double t66;
  double t69;
  double t7;
  double t70;
  double t71;
  double t73;
  double t74;
  double t79;
  double t8;
  double t82;
  double t84;
  double t86;
  double t88;
  double t9;
  double t90;
  double t92;
  double t99;
  t2 = 0.1e1 / (0.4e1 - s);
  t3 = s * t2;
  t4 = sp(s1, k2);
  t7 = beta_y * s;
  t8 = sp(s1, p1);
  t9 = t2 * t8;
  t10 = t7 * t9;
  t12 = -0.64e2 * t3 * t4 - 0.32e2 * t10;
  t13 = sp(k1, s2);
  t15 = t2 * t4;
  t16 = sp(s2, p1);
  t18 = t7 * t15 * t16;
  t19 = 0.32e2 * t18;
  t20 = beta * beta;
  t22 = y * y;
  t23 = t20 * s * t22;
  t26 = s * (t23 - 0.2e1 * s + 0.8e1);
  t27 = sp(s1, s2);
  t28 = t2 * t27;
  t30 = 0.16e2 * t26 * t28;
  t32 = 0.16e2 * t26 * t2;
  t34 = RE(I2_MT2_0_MT2_MU2_0);
  t38 = -t12 * t13 - t19 + t30 - t32;
  t39 = RE(I2_S12_0_0_MU2_0);
  t41 = 0.3e1 * s;
  t43 = s * (0.4e1 - t41);
  t51 = s * s;
  t55 = s * (0.4e1 * t23 - 0.3e1 * t51 + 0.16e2 * s - 0.16e2);
  t60 = (-0.16e2 * t15 * t43 + 0.64e2 * t10) * t13 - 0.64e2 * t18 + 0.8e1 * t55 * t28 - 0.8e1 * t55 * t2;
  t61 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t63 = EPS_(k1, k2, p1, s1);
  t66 = EPS_(k1, k2, p1, s2);
  t69 = -t2 * t63 * t7 - t2 * t66 * t7;
  t70 = 0.32e2 * t69;
  t71 = IM(I2_S12_0_0_MU2_0);
  t73 = 0.64e2 * t69;
  t74 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t79 = EPS_(k1, k2, s1, s2);
  t82 = EPS_(k1, p1, s1, s2);
  t84 = t7 * t2 * t82;
  t86 = EPS_(k2, p1, s1, s2);
  t88 = t7 * t2 * t86;
  t90 = -0.64e2 * t3 * t79 - 0.32e2 * t84 - 0.32e2 * t88;
  t92 = -t90;
  t99 = -0.16e2 * t2 * t43 * t79 + 0.64e2 * t84 + 0.64e2 * t88;
  t103 = t2 * t13;
  t110 = beta_y * t51;
  t111 = t110 * t9;
  t114 = t110 * t2 * t16;
  t116 = 0.8e1 * s * (t23 + t7 + 0.8e1) * t103 + 0.8e1 * s * (t23 - t7 + 0.8e1) * t15 + 0.16e2 * t111 - 0.16e2 * t114;
  t128 = 0.16e2 * s * (t23 + t7 + t41 - 0.4e1) * t103 + 0.16e2 * s * (t23 - t7 + t41 - 0.4e1) * t15 + 0.32e2 * t111 - 0.32e2 * t114;
  return(At_fH_re * ((t12 * t13 + t19 - t30 + t32) * t34 + t38 * t39 + t60 * t61 + t70 * t71 + t73 * t74) * PREF_V_CA + Bt_fH_re * (t116 * t71 + t128 * t74 + t34 * t90 + t39 * t92 + t61 * t99) * PREF_V_CA + At_fH_im * (t34 * t70 + t38 * t71 - t39 * t70 + t60 * t74 - t61 * t73) * PREF_V_CA + Bt_fH_im * (t116 * t34 - t116 * t39 - t128 * t61 + t71 * t92 + t74 * t99) * PREF_V_CA);
}