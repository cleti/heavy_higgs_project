
#include "../../../../inc/Functions_Shared.h"

double Eval_V_B3 (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  c_double& cg  = I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_0;
  c_double& cg1 = I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_1;
  // c_double& cg3 = I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_0;
  // c_double& cg5 = I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_1;

  double t1;
  double t10;
  double t102;
  double t105;
  double t108;
  double t111;
  double t113;
  double t116;
  double t117;
  double t118;
  double t12;
  double t13;
  double t139;
  double t14;
  double t141;
  double t143;
  double t15;
  double t153;
  double t154;
  double t156;
  double t158;
  double t159;
  double t16;
  double t161;
  double t162;
  double t163;
  double t164;
  double t165;
  double t166;
  double t171;
  double t18;
  double t184;
  double t196;
  double t22;
  double t25;
  double t28;
  double t29;
  double t30;
  double t31;
  double t32;
  double t38;
  double t4;
  double t41;
  double t42;
  double t45;
  double t47;
  double t48;
  double t5;
  double t50;
  double t52;
  double t54;
  double t55;
  double t56;
  double t59;
  double t6;
  double t60;
  double t61;
  double t64;
  double t69;
  double t7;
  double t73;
  double t77;
  double t78;
  double t81;
  c_double t82;
  double t83;
  double t85;
  double t88;
  double t93;
  double t99;
  t1 = beta_y * beta_y;
  t4 = s * s;
  t5 = t4 * s;
  t6 = (0.3e1 * t1 + 0.5e1) * t5;
  t7 = RE(cg);
  t10 = beta_y + 0.1e1;
  t12 = beta_y - 0.1e1;
  t13 = t12 * t12;
  t14 = t1 * s;
  t15 = t14 - s + 0.4e1;
  t16 = 0.1e1 / t15;
  t18 = RE(I3_0_T11_MT2_0_0_MT2_MU2_1);
  t22 = 0.3e1 * beta_y;
  t25 = RE(I3_0_T12_MT2_0_0_MT2_MU2_0);
  t28 = t4 * t4;
  t29 = t28 * t13;
  t30 = t10 * t10;
  t31 = t30 * t16;
  t32 = RE(cg1);
  t38 = RE(I3_0_T11_MT2_0_0_MT2_MU2_0);
  t41 = t4 * t10;
  t42 = RE(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
  t45 = t1 * t1;
  t47 = 0.2e1 * t45 * s;
  t48 = 0.5e1 * t14;
  t50 = 0.3e1 * s;
  t52 = 0.1e1 / t30;
  t54 = 0.1e1 / t13;
  t55 = RE(I2_MT2_0_MT2_MU2_0);
  t56 = t54 * t55;
  t59 = beta_y * s;
  t60 = 0.3e1 * t59;
  t61 = 0.2e1 * s;
  t64 = RE(I2_T11_0_MT2_MU2_0);
  t69 = RE(I2_T12_0_MT2_MU2_0);
  t73 = RE(I3_0_T12_MT2_0_0_MT2_MU2_1);
  t77 = t12 * t4;
  t78 = RE(I3_MT2_0_T11_0_MT2_MT2_MU2_0);
  t81 = 0.2e1 * t6 * t7 + 0.12e2 * t5 * t10 * t13 * t16 * t18 - 0.4e1 * (t22 - 0.1e1) * t4 * t25 + 0.6e1 * t29 * t31 * t32 + 0.4e1 * (t22 + 0.1e1) * t4 * t38 - 0.12e2 * t41 * t42 - 0.32e2 * (t47 - t48 - 0.4e1 * t1 + t50 - 0.4e1) * t52 * t56 + 0.32e2 * (t14 - t60 + t61 - 0.2e1) * t54 * t64 + 0.32e2 * (t14 + t60 + t61 - 0.2e1) * t52 * t69 - 0.12e2 * t5 * t12 * t31 * t73 + 0.12e2 * t77 * t78;
  t82 = DenS(s, mH, GammaH);
  t83 = RE(t82);
  t85 = IM(cg);
  t88 = IM(cg1);
  t93 = IM(t82);
  t99 = AlphaS3 / 0.3141592653589793e1;
  t102 = VF(FA0 * CF * CA2 * t99 * Bt);
  t105 = 0.8e1 * t1;
  t108 = (-0.3e1 * t14 + t105 - 0.5e1 * s + 0.24e2) * t4;
  t111 = t5 * t13 * t30;
  t113 = (0.8e1 - t50) * t16;
  t116 = 0.2e1 * t14;
  t117 = 0.7e1 * t59;
  t118 = 0.8e1 * beta_y;
  t139 = 0.1e1 / (t59 - s + 0.2e1);
  t141 = 0.1e1 / (t59 + s - 0.2e1);
  t143 = RE(I1_MT2_MU2_0);
  t153 = beta_y * t1;
  t154 = t153 * t4;
  t156 = 0.5e1 * t153 * s;
  t158 = 0.4e1 * t1 * t4;
  t159 = 0.17e2 * t14;
  t161 = 0.5e1 * beta_y * t4;
  t162 = 0.23e2 * t59;
  t163 = 0.2e1 * t4;
  t164 = 0.16e2 * beta_y;
  t165 = 0.11e2 * s;
  t166 = -t154 + t156 + t158 - t159 - t161 + t105 + t162 + t163 - t164 - t165 + 0.12e2;
  t171 = -t154 + t156 - t158 + t159 - t161 - t105 + t162 - t163 - t164 + t165 - 0.12e2;
  t184 = t108 * t7 + t111 * t113 * t32 + 0.2e1 * (t116 + t117 - t118 + s) * s * t25 - 0.2e1 * t77 * t30 * t113 * t73 + 0.2e1 * (t116 - t117 + t118 + s) * s * t38 + 0.2e1 * t41 * t13 * t113 * t18 - 0.16e2 * (0.2e1 - s) * t15 * s * t139 * t141 * t143 - 0.16e2 * (-t47 + 0.8e1 * t45 + t48 - 0.12e2 * t1 - t50 + 0.12e2) * t52 * t56 + 0.16e2 * t166 * t139 * t54 * t64 + 0.16e2 * t171 * t52 * t141 * t69 + 0.2e1 * (-t60 + t118 + t50 - 0.16e2) * s * t78 - 0.2e1 * (-t60 + t118 - t50 + 0.16e2) * s * t42;
  t196 = VF(FH0 * CF * CA2 * t99 * At);
  return((t81 * t83 + (0.6e1 * t29 * t31 * t88 + 0.2e1 * t6 * t85) * t93) * t102 + (t184 * t83 + (t111 * t113 * t88 + t108 * t85) * t93) * t196);
}
