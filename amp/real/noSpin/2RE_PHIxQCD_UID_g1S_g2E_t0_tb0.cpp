
#include "UID_HEADER.h"

double Eval_UID_SE00 (UID_ARGS)
{

UID_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);


  double t10;
  double t102;
  double t113;
  double t118;
  double t119;
  double t122;
  double t13;
  double t147;
  double t15;
  double t150;
  double t18;
  double t19;
  double t20;
  double t24;
  double t25;
  double t28;
  double t3;
  double t30;
  double t38;
  double t39;
  double t4;
  double t40;
  double t44;
  double t46;
  double t47;
  double t48;
  double t5;
  double t50;
  double t51;
  double t53;
  double t54;
  double t57;
  double t62;
  double t64;
  double t66;
  double t69;
  double t7;
  double t71;
  double t73;
  double t76;
  double t8;
  double t86;
  double t88;
  double t95;
  double t96;
  double t98;
  double t99;
  t3 = sp(p1, p2);
  t4 = 0.1e1 / t3;
  t5 = sp(p1, p3);
  t7 = sp(p2, p3);
  t8 = t4 * t5 * t7;
  t10 = sp(K2, p3);
  t13 = t5 * t5;
  t15 = sp(K1, p3);
  t18 = sp(K2, p1);
  t19 = t4 * t18;
  t20 = t5 * t7;
  t24 = sp(K1, p1);
  t25 = 0.1e1 / t24;
  t28 = sp(K1, P2);
  t30 = sp(P2, p3);
  t38 = sp(P2, p1);
  t39 = t4 * t38;
  t40 = t39 * t7;
  t44 = t30 * t7;
  t46 = 0.512e3 * t19 * t44;
  t47 = sp(K2, P2);
  t48 = t47 * t13;
  t50 = t18 * t30;
  t51 = t50 * t5;
  t53 = t10 * t38;
  t54 = 0.512e3 * t53;
  t57 = t38 * t7;
  t62 = t53 * t5;
  t64 = t15 * t15;
  t66 = 0.512e3 * t64 * t38;
  t69 = t30 * t5;
  t71 = t39 * t20;
  t73 = -t71 - t69;
  t76 = t24 * t24;
  t86 = t39 * t44;
  t88 = t30 * t30;
  t95 = t38 * t30;
  t96 = t38 * t38;
  t98 = t4 * t96 * t7;
  t99 = -t95 + t98;
  t102 = t88 - t86;
  t113 = 0.1e1 / t28;
  t118 = EPS_(K1, K2, p1, p3);
  t119 = t118 * t4;
  t122 = t118 * t15;
  t147 = -0.512e3 * t99;
  t150 = -0.512e3 * t102;
  return(At_fH_re * ((-0.512e3 * t8 + (-0.512e3 * t10 * t5 + 0.512e3 * t15 * t5 + 0.512e3 * t19 * t20 - 0.256e3 * t13) * t25) * t28 - 0.512e3 * t10 * t30 + (-0.512e3 * t30 * t4 * t7 - 0.512e3 * t8) * t24 + (0.512e3 * t40 + 0.512e3 * t30 + 0.512e3 * t5) * t15 + t46 + (0.256e3 * t48 - 0.256e3 * t51 + (-0.512e3 * t19 * t57 + 0.256e3 * t38 * t5 + t54) * t15 - 0.256e3 * t62 - t66) * t25 + 0.512e3 * t69 + 0.256e3 * t71 + (-t66 + 0.256e3 * t73 * t47 - 0.512e3 * t4 * t76 * t44 + (0.512e3 * t4 * t47 * t20 - 0.512e3 * t4 * t10 * t57 + 0.256e3 * t86 + t46 - 0.256e3 * t88 + (0.512e3 * t40 + 0.512e3 * t30) * t15) * t24 + 0.256e3 * t99 * t10 + 0.256e3 * t102 * t18 + (-0.512e3 * t47 * t5 - 0.512e3 * t50 + t54 + 0.256e3 * t95 - 0.256e3 * t98) * t15) * t113) * PREF_R_CA + At_fA_re * (-0.1024e4 * t119 * t7 + (-0.2048e4 * t119 * t24 * t7 - 0.1024e4 * t118 * t30 + 0.1024e4 * t119 * t57 + 0.2048e4 * t122) * t113 + (0.1024e4 * t119 * t18 * t7 - 0.1024e4 * t10 * t118 - 0.1024e4 * t118 * t5 + 0.1024e4 * t122) * t25) * PREF_R_CA + Bt_fA_re * (-0.512e3 * t25 * t28 * t13 + (t10 * t147 + t147 * t15 + t150 * t18 + t150 * t24 - 0.512e3 * t73 * t47) * t113 + 0.512e3 * t71 + (0.512e3 * t15 * t38 * t5 - 0.512e3 * t48 + 0.512e3 * t51 + 0.512e3 * t62) * t25 + 0.1024e4 * t69) * PREF_R_CA);
}
