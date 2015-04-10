
#include "UID_HEADER.h"

double Eval_UID_PHIxPHI_ES00 (UID_ARGS)
{

UID_DEFINITIONS


  double t101;
  double t104;
  double t11;
  double t12;
  double t120;
  double t125;
  double t127;
  double t13;
  double t134;
  double t14;
  double t15;
  double t16;
  double t18;
  double t19;
  double t27;
  double t28;
  double t3;
  double t37;
  double t4;
  double t44;
  double t45;
  double t47;
  double t48;
  double t50;
  double t51;
  double t53;
  double t54;
  double t57;
  double t59;
  double t6;
  double t60;
  double t62;
  double t63;
  double t66;
  double t68;
  double t69;
  double t7;
  double t70;
  double t72;
  double t73;
  double t76;
  double t80;
  double t82;
  double t84;
  double t85;
  double t87;
  double t9;
  double t99;
  t3 = sp(p2, p3);
  t4 = sp(P1, p3);
  t6 = sp(P1, p2);
  t7 = sp(K2, S1);
  t9 = sp(K1, S2);
  t11 = t3 * t4 * t6 * t7 * t9;
  t12 = 0.1024e4 * t11;
  t13 = sp(S1, S2);
  t14 = t3 * t13;
  t15 = -t14 + t3;
  t16 = 0.1024e4 * t15;
  t18 = sp(K1, K2);
  t19 = t6 * t18;
  t27 = 0.4096e4 * t11;
  t28 = 0.4096e4 * t15;
  t37 = EPS_(K1, K2, S1, S2);
  t44 = EPS_(K1, S2, p2, p3);
  t45 = EPS_(K2, S2, p2, p3);
  t47 = -0.4096e4 * t44 - 0.4096e4 * t45;
  t48 = sp(S1, p3);
  t50 = EPS_(K1, S1, p2, p3);
  t51 = EPS_(K2, S1, p2, p3);
  t53 = 0.4096e4 * t50 + 0.4096e4 * t51;
  t54 = sp(S2, p3);
  t57 = sp(K2, P1);
  t59 = -t47;
  t60 = sp(S1, P1);
  t62 = -t53;
  t63 = sp(S2, P1);
  t66 = sp(K2, p3);
  t68 = EPS_(K1, K2, p2, p3);
  t69 = t68 * t54;
  t70 = t69 * t60;
  t72 = t68 * t48;
  t73 = t72 * t63;
  t76 = sp(K1, p2);
  t80 = sp(K2, p2);
  t82 = sp(S1, p2);
  t84 = t68 * t82 * t63;
  t85 = sp(S2, p2);
  t87 = t68 * t85 * t60;
  t99 = t69 * t82;
  t101 = t72 * t85;
  t104 = sp(K1, P1);
  t120 = sp(K1, p3);
  t125 = t14 + t3;
  t127 = 0.1024e4 * t125 * t4;
  t134 = 0.4096e4 * t125 * t4;
  return(At2_fH2_De * (t16 * t19 * t4 - t16 * t4 * t6 + t12) * PREF_R_PHI_CA + At2_fA2_De * (t19 * t28 * t4 - t28 * t4 * t6 + t27) * PREF_R_PHI_CA + 0.2048e4 * At_Bt_fH2_De * t37 * t3 * t4 * t6 * PREF_R_PHI_CA + At_Bt_fA2_De * (((t47 * t48 + t53 * t54) * t57 + (t59 * t60 + t62 * t63) * t66 - 0.4096e4 * t70 + 0.4096e4 * t73) * t76 + (-0.4096e4 * t70 + 0.4096e4 * t73) * t80 + (-0.4096e4 * t84 + 0.4096e4 * t87) * t66 + ((t48 * t59 + t54 * t62) * t80 + (t47 * t82 + t53 * t85) * t66 + 0.4096e4 * t99 - 0.4096e4 * t101) * t104 + (0.4096e4 * t99 - 0.4096e4 * t101) * t57 + ((t59 * t82 + t62 * t85) * t57 + (t47 * t60 + t53 * t63) * t80 + 0.4096e4 * t87 - 0.4096e4 * t84) * t120) * PREF_R_PHI_CA + Bt2_fH2_De * (t127 * t19 + t127 * t6 - t12) * PREF_R_PHI_CA + Bt2_fA2_De * (t134 * t19 + t134 * t6 - t27) * PREF_R_PHI_CA);
}
