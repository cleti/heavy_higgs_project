
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR (AMP_ARGS)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);


  double t100;
  double t102;
  double t103;
  double t104;
  double t106;
  double t108;
  double t11;
  double t116;
  double t119;
  double t130;
  double t132;
  double t133;
  double t14;
  double t16;
  double t22;
  double t24;
  double t25;
  double t27;
  double t28;
  double t29;
  double t3;
  double t32;
  double t34;
  double t35;
  double t38;
  double t4;
  double t41;
  double t45;
  double t47;
  double t48;
  double t49;
  double t5;
  double t50;
  double t51;
  double t52;
  double t55;
  double t56;
  double t57;
  double t58;
  double t59;
  double t6;
  double t63;
  double t66;
  double t67;
  double t68;
  double t69;
  double t7;
  double t70;
  double t71;
  double t73;
  double t74;
  double t75;
  double t8;
  double t80;
  double t83;
  double t85;
  double t87;
  double t89;
  double t9;
  double t98;
  double t99;
  t3 = -0.2048e4 * At_Bt_fA2_De - 0.512e3 * At_Bt_fH2_De;
  t4 = t3 * PREF_R_PHI_CF;
  t5 = sp(p1, p2);
  t6 = t5 * t5;
  t7 = sp(k2, p3);
  t8 = t7 * t7;
  t9 = 0.1e1 / t8;
  t11 = t4 * t6 * t9;
  t14 = 0.4096e4 * At_Bt_fA2_De + 0.1024e4 * At_Bt_fH2_De;
  t16 = t6 * t5;
  t22 = 0.1e1 / t7;
  t24 = sp(k1, p3);
  t25 = 0.1e1 / t24;
  t27 = t24 * t24;
  t28 = 0.1e1 / t27;
  t29 = t6 * t28;
  t32 = EPS_(k1, k2, s1, s2);
  t34 = t6 * t22;
  t35 = t34 * t25;
  t38 = EPS_(k1, p3, s1, s2);
  t41 = -t3 * PREF_R_PHI_CF;
  t45 = EPS_(k2, p3, s1, s2);
  t47 = 0.1024e4 * At2_fA2_De;
  t48 = 0.256e3 * At2_fH2_De;
  t49 = 0.1024e4 * Bt2_fA2_De;
  t50 = 0.256e3 * Bt2_fH2_De;
  t51 = t47 + t48 + t49 + t50;
  t52 = t51 * PREF_R_PHI_CF;
  t55 = -t47 - t48 + t49 + t50;
  t56 = t55 * PREF_R_PHI_CF;
  t57 = sp(s1, k2);
  t58 = t56 * t57;
  t59 = sp(s1, p3);
  t63 = sp(k1, s2);
  t66 = 0.2048e4 * At2_fA2_De;
  t67 = 0.512e3 * At2_fH2_De;
  t68 = 0.2048e4 * Bt2_fA2_De;
  t69 = 0.512e3 * Bt2_fH2_De;
  t70 = t66 + t67 + t68 + t69;
  t71 = t70 * PREF_R_PHI_CF;
  t73 = t66 + t67 - t68 - t69;
  t74 = t73 * PREF_R_PHI_CF;
  t75 = sp(s1, s2);
  t80 = (-t70 * PREF_R_PHI_CF + t74 * t75) * t16;
  t83 = -0.4096e4 * At2_fA2_De - 0.1024e4 * At2_fH2_De;
  t85 = t83 * PREF_R_PHI_CF * t75;
  t87 = -t83 * PREF_R_PHI_CF;
  t89 = (t85 + t87) * t6;
  t98 = (-t55 * t75 * PREF_R_PHI_CF - t51 * PREF_R_PHI_CF) * t16;
  t99 = sp(s2, p3);
  t100 = t99 * t57;
  t102 = t99 * t59;
  t103 = t52 * t102;
  t104 = -t66 - t67;
  t106 = t104 * PREF_R_PHI_CF * t75;
  t108 = -t104 * PREF_R_PHI_CF;
  t116 = -t73 * PREF_R_PHI_CF;
  t119 = (t68 + t69) * PREF_R_PHI_CF;
  t130 = t6 * t6;
  t132 = 0.6144e4 * At2_fA2_De;
  t133 = 0.1536e4 * At2_fH2_De;
  return((t11 + (t14 * t16 * PREF_R_PHI_CF - t14 * t6 * PREF_R_PHI_CF) * t22 * t25 + t4 * t29) * t32 + (t35 * t4 + t11) * t38 + (t29 * t41 + t35 * t41) * t45 + t52 * t34 * t24 + (t56 * t59 + t58) * t6 * t9 * t63 + t71 * t6 + (t80 + t89) * t22 + (t98 + (t100 * t52 + t103 + t106 + t108) * t6) * t9 + ((t74 * t57 * t16 + (t116 * t57 + t119 * t59) * t6) * t22 * t63 + t52 * t6 * t7 + t80 + t89 + ((t116 * t75 + t71) * t130 + ((t132 + t133 - t68 - t69) * PREF_R_PHI_CF * t75 + (-t132 - t133 - t68 - t69) * PREF_R_PHI_CF) * t16 + (t100 * t119 + t102 * t71 + t85 + t87) * t6) * t22) * t25 + ((t52 * t59 + t58) * t6 * t63 + t98 + (t100 * t56 + t103 + t106 + t108) * t6) * t28);
}
