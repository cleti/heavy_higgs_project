
#include "AMP_HEADER.h"

double Eval_V_B1 (AMP_ARGS)
{

  AMP_DEFINITIONS
  HP_REFS_PHIxQCD(hp)
  AP_REFS_V(ap)

    
  c_double& cg  = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg1 = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1;
  c_double& cg3 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg5 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1;

  
  double t1;
  double t106;
  double t109;
  double t11;
  double t112;
  double t113;
  double t114;
  double t12;
  double t120;
  double t127;
  double t128;
  double t13;
  double t131;
  double t134;
  double t138;
  double t142;
  double t146;
  double t18;
  double t21;
  double t24;
  double t25;
  double t27;
  double t28;
  double t29;
  double t3;
  double t30;
  double t31;
  double t32;
  double t36;
  double t37;
  double t4;
  double t44;
  double t45;
  double t49;
  double t50;
  double t54;
  double t58;
  double t61;
  double t64;
  double t65;
  double t66;
  double t67;
  double t68;
  double t69;
  double t72;
  double t75;
  double t76;
  double t8;
  double t80;
  double t81;
  double t82;
  double t83;
  double t87;
  double t9;
  double t91;
  double t92;
  double t93;
  double t97;
  t1 = beta * beta;
  t3 = y * y;
  t4 = t1 * s * t3;
  t8 = beta_y + 0.1e1;
  t9 = 0.1e1 / t8;
  t11 = beta_y - 0.1e1;
  t12 = 0.1e1 / t11;
  t13 = RE(I2_MT2_0_MT2_MU2_0);
  t18 = s * s;
  t21 = t1 * beta * t18 * t3 * y;
  t24 = 0.3e1 * t1 * t18 * t3;
  t25 = 0.8e1 * t4;
  t27 = 0.9e1 * beta_y * t18;
  t28 = beta_y * s;
  t29 = 0.36e2 * t28;
  t30 = 0.5e1 * t18;
  t31 = 0.16e2 * beta_y;
  t32 = 0.28e2 * s;
  t36 = 0.1e1 / (t28 - s + 0.2e1);
  t37 = RE(I2_T11_0_MT2_MU2_0);
  t44 = 0.1e1 / (t28 + s - 0.2e1);
  t45 = RE(I2_T12_0_MT2_MU2_0);
  t49 = 0.2e1 - s;
  t50 = t4 - s + 0.4e1;
  t54 = RE(I1_MT2_MU2_0);
  t58 = RE(I3_MT2_0_T11_0_MT2_MT2_MU2_0);
  t61 = RE(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
  t64 = -s + 0.4e1;
  t65 = t64 * t64;
  t66 = s * t65;
  t67 = 0.1e1 / t50;
  t68 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t69 = t67 * t68;
  t72 = RE(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t75 = t49 * t64;
  t76 = RE(cg);
  t80 = t18 * t65;
  t81 = t11 * t11;
  t82 = t81 * t67;
  t83 = RE(cg1);
  t87 = RE(cg3);
  t91 = t8 * t8;
  t92 = t91 * t67;
  t93 = RE(cg5);
  t97 = 0.64e2 * (0.2e1 * t1 * t3 + s - t4 - 0.4e1) * t9 * t12 * t13 + 0.8e1 * (t21 + t24 - t25 - t27 + t29 + t30 - t31 - t32 + 0.32e2) * t12 * t36 * t37 - 0.8e1 * (t21 - t24 + t25 - t27 + t29 - t30 - t31 + t32 - 0.32e2) * t9 * t44 * t45 + 0.16e2 * t49 * t50 * s * t36 * t44 * t54 + 0.32e2 * s * t58 + 0.32e2 * s * t61 - 0.32e2 * t66 * t69 - 0.64e2 * s * t72 - 0.8e1 * t75 * s * t76 - 0.8e1 * t80 * t82 * t83 - 0.8e1 * t75 * s * t87 - 0.8e1 * t80 * t92 * t93;
  t106 = t18 * t64;
  t109 = t49 * t18;
  t112 = t18 * s;
  t113 = t112 * t81;
  t114 = t64 * t67;
  t120 = t112 * t91;
  t127 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t128 = t67 * t127;
  t131 = IM(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t134 = IM(cg);
  t138 = IM(cg1);
  t142 = IM(cg3);
  t146 = IM(cg5);
  return(At_fH_re * t97 * PREF_V_CF + Bt_fA_re * (-0.16e2 * t113 * t114 * t83 - 0.16e2 * t114 * t120 * t93 + 0.64e2 * s * t13 - 0.32e2 * s * t37 - 0.32e2 * s * t45 - 0.64e2 * t106 * t69 - 0.16e2 * t109 * t76 - 0.16e2 * t109 * t87) * PREF_V_CF + At_fH_im * (-0.8e1 * s * t134 * t75 - 0.8e1 * s * t142 * t75 - 0.8e1 * t138 * t80 * t82 - 0.8e1 * t146 * t80 * t92 - 0.64e2 * s * t131 - 0.32e2 * t128 * t66) * PREF_V_CF + Bt_fA_im * (-0.16e2 * t113 * t114 * t138 - 0.16e2 * t114 * t120 * t146 - 0.64e2 * t106 * t128 - 0.16e2 * t109 * t134 - 0.16e2 * t109 * t142) * PREF_V_CF);
}
