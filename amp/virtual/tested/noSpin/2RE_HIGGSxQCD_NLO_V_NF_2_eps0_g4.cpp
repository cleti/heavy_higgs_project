
#include "../../../../inc/Functions_Shared.h"

double Eval_V_NF_2 (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  double t1;
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t21;
  double t26;
  double t27;
  double t28;
  c_double t3;
  double t30;
  double t31;
  double t35;
  double t4;
  double t40;
  double t42;
  c_double t43;
  double t44;
  double t46;
  double t49;
  double t5;
  double t50;
  double t55;
  double t57;
  double t60;
  double t62;
  double t67;
  double t68;
  double t7;
  double t71;
  double t74;
  double t77;
  double t79;
  double t8;
  double t82;
  double t84;
  double t9;
  double t90;
  double t92;
  double t98;
  t1 = s * s;
  t2 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_0);
  t3 = MH2;
  t4 = RE(t3);
  t5 = t2 * t4;
  t7 = 0.2e1 * t5 * t1;
  t8 = IM(I2_S12_0_MH2_MU2_0);
  t9 = IM(t3);
  t11 = t8 * t9 * s;
  t12 = RE(I2_MT2_MT2_MH2_MU2_0);
  t13 = t12 * t4;
  t14 = t13 * s;
  t15 = IM(I2_MT2_MT2_MH2_MU2_0);
  t16 = t15 * t9;
  t17 = t16 * s;
  t18 = RE(I2_S12_0_MH2_MU2_0);
  t20 = t18 * t4 * s;
  t21 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_2);
  t26 = 0.4e1 * t5 * s;
  t27 = IM(I3_MT2_MT2_S12_0_MT2_MH2_MU2_0);
  t28 = t27 * t9;
  t30 = 0.4e1 * t28 * s;
  t31 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_1);
  t35 = IM(I3_MT2_MT2_S12_0_MT2_MH2_MU2_1);
  t40 = 0.2e1 * t28 * t1;
  t42 = 0.2e1 * t12 * s;
  t43 = t3 * t3;
  t44 = RE(t43);
  t46 = 0.2e1 * t2 * t44;
  t49 = 0.4e1 * s * t21 * t4 + 0.4e1 * s * t31 * t4 - 0.4e1 * s * t35 * t9 - 0.2e1 * t31 * t44 - t1 - t11 - t14 + t17 + t20 - t26 + t30 - t40 + t42 - t46 + t7;
  t50 = IM(t43);
  t55 = RE(I2_MT2_0_MT2_MU2_0);
  t57 = 0.2e1 * t55 * t4;
  t60 = 0.2e1 * t13;
  t62 = 0.6e1 * t55 * s;
  t67 = 0.2e1 * t16;
  t68 = t2 * t1;
  t71 = 0.2e1 * t27 * t50;
  t74 = t18 * t1;
  t77 = 0.2e1 * t55 * t1;
  t79 = 0.8e1 * s * t18 + 0.32e2 * s * t2 - 0.2e1 * t1 * t21 - 0.2e1 * t1 * t31 - 0.2e1 * t21 * t44 + 0.2e1 * t35 * t50 + 0.4e1 * s - t57 + t60 - t62 - t67 - 0.10e2 * t68 + t71 - 0.3e1 * t74 + t77;
  t82 = beta_y * beta_y;
  t84 = t82 / s;
  t90 = 0.1e1 / (0.4e1 - s) / (beta_y - 0.1e1);
  t92 = 0.1e1 / (beta_y + 0.1e1);
  t98 = -t26 - 0.2e1 * t68 + t30 + t7 + t60 - t57 - t46 - t40 + t71 - t67 + t42 - t62 + t20 - t74 - t14 - t11 + t17 + t77;
  return(0.16e2 * FH0 * (t49 + t79) * t84 * t90 * t92 * PREF_V_CA_At - 0.32e2 * FA0 * t98 * t84 * t90 * t92 * PREF_V_CA_Bt);
}
