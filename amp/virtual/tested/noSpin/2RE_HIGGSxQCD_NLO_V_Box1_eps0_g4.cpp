
#include "../../../../inc/Functions_Shared.h"

double Eval_V_B1 (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  c_double& cg  = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg1 = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1;
  c_double& cg3 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg5 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1;


  double t1;
  double t10;
  double t100;
  double t102;
  double t103;
  double t105;
  double t106;
  double t107;
  double t108;
  double t109;
  double t11;
  double t110;
  double t114;
  double t12;
  double t121;
  double t125;
  double t128;
  double t13;
  double t134;
  double t138;
  double t14;
  double t142;
  double t149;
  double t15;
  double t156;
  double t16;
  double t17;
  double t175;
  double t2;
  double t21;
  double t22;
  double t23;
  double t24;
  double t28;
  double t29;
  double t3;
  double t30;
  double t33;
  double t36;
  double t39;
  double t4;
  double t42;
  c_double t46;
  double t47;
  double t49;
  double t52;
  double t56;
  double t59;
  double t63;
  double t64;
  double t68;
  double t7;
  double t74;
  double t77;
  double t79;
  double t8;
  double t82;
  double t83;
  double t84;
  double t88;
  double t9;
  double t93;
  double t95;
  t1 = 0.2e1 - s;
  t2 = s * s;
  t3 = t1 * t2;
  t4 = RE(cg);
  t7 = t2 * s;
  t8 = beta_y - 0.1e1;
  t9 = t8 * t8;
  t10 = t7 * t9;
  t11 = 0.4e1 - s;
  t12 = beta_y * beta_y;
  t13 = t12 * s;
  t14 = t13 - s + 0.4e1;
  t15 = 0.1e1 / t14;
  t16 = t11 * t15;
  t17 = RE(cg1);
  t21 = beta_y + 0.1e1;
  t22 = t21 * t21;
  t23 = t7 * t22;
  t24 = RE(cg5);
  t28 = t2 * t11;
  t29 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t30 = t15 * t29;
  t33 = RE(I2_MT2_0_MT2_MU2_0);
  t36 = RE(I2_T11_0_MT2_MU2_0);
  t39 = RE(I2_T12_0_MT2_MU2_0);
  t42 = RE(cg3);
  t46 = DenS(s, mH, GammaH);
  t47 = RE(t46);
  t49 = IM(cg3);
  t52 = IM(cg5);
  t56 = IM(cg);
  t59 = IM(cg1);
  t63 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t64 = t15 * t63;
  t68 = IM(t46);
  t74 = AlphaS3 / 0.3141592653589793e1;
  t77 = VF(FA0 * CF2 * CA * t74 * Bt);
  t79 = RE(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
  t82 = t11 * t11;
  t83 = t2 * t82;
  t84 = t22 * t15;
  t88 = RE(I3_MT2_0_T11_0_MT2_MT2_MU2_0);
  t93 = 0.1e1 / t8;
  t95 = 0.1e1 / t21;
  t100 = t12 * beta_y * t2;
  t102 = 0.3e1 * t12 * t2;
  t103 = 0.8e1 * t13;
  t105 = 0.9e1 * beta_y * t2;
  t106 = beta_y * s;
  t107 = 0.36e2 * t106;
  t108 = 0.5e1 * t2;
  t109 = 0.16e2 * beta_y;
  t110 = 0.28e2 * s;
  t114 = 0.1e1 / (t106 - s + 0.2e1);
  t121 = 0.1e1 / (t106 + s - 0.2e1);
  t125 = s * t82;
  t128 = RE(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t134 = RE(I1_MT2_MU2_0);
  t138 = t1 * t11;
  t142 = t9 * t15;
  t149 = 0.32e2 * s * t79 - 0.8e1 * t83 * t84 * t24 + 0.32e2 * s * t88 + 0.64e2 * (-t13 + 0.2e1 * t12 + s - 0.4e1) * t93 * t95 * t33 + 0.8e1 * (t100 + t102 - t103 - t105 + t107 + t108 - t109 - t110 + 0.32e2) * t93 * t114 * t36 - 0.8e1 * (t100 - t102 + t103 - t105 + t107 - t108 - t109 + t110 - 0.32e2) * t95 * t121 * t39 - 0.32e2 * t125 * t30 - 0.64e2 * s * t128 + 0.16e2 * t1 * t14 * s * t121 * t114 * t134 - 0.8e1 * t138 * s * t4 - 0.8e1 * t83 * t142 * t17 - 0.8e1 * t138 * s * t42;
  t156 = IM(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t175 = VF(FH0 * CF2 * CA * t74 * At);
  return(((-0.16e2 * t10 * t16 * t17 - 0.16e2 * t16 * t23 * t24 + 0.64e2 * s * t33 - 0.32e2 * t36 * s - 0.32e2 * s * t39 - 0.64e2 * t28 * t30 - 0.16e2 * t3 * t4 - 0.16e2 * t3 * t42) * t47 + (-0.16e2 * t10 * t16 * t59 - 0.16e2 * t16 * t23 * t52 - 0.64e2 * t28 * t64 - 0.16e2 * t3 * t49 - 0.16e2 * t3 * t56) * t68) * t77 + (t149 * t47 + (-0.8e1 * s * t138 * t49 - 0.8e1 * s * t138 * t56 - 0.8e1 * t142 * t59 * t83 - 0.8e1 * t52 * t83 * t84 - 0.64e2 * s * t156 - 0.32e2 * t125 * t64) * t68) * t175);
}