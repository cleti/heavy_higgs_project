
#include "../../../../inc/Functions_Shared.h"

double Eval_V_NF_3 (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  c_double t1;
  double t10;
  double t12;
  double t14;
  double t15;
  double t17;
  double t19;
  double t2;
  double t22;
  double t23;
  double t25;
  c_double t26;
  double t27;
  double t29;
  double t3;
  double t30;
  double t32;
  double t33;
  double t34;
  double t35;
  double t36;
  double t39;
  double t4;
  double t41;
  double t42;
  double t43;
  double t44;
  double t46;
  double t47;
  double t48;
  double t50;
  double t52;
  double t58;
  double t60;
  double t65;
  double t73;
  double t77;
  double t8;
  double t81;
  double t9;
  double t94;
  t1 = MH2;
  t2 = RE(t1);
  t3 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_0);
  t4 = t2 * t3;
  t8 = IM(t1);
  t9 = IM(I3_MT2_MT2_S12_0_MT2_MH2_MU2_0);
  t10 = t8 * t9;
  t12 = s * s;
  t14 = 0.2e1 * t3 * t12;
  t15 = RE(I2_MT2_MT2_MH2_MU2_0);
  t17 = RE(I2_MT2_0_MT2_MU2_0);
  t19 = RE(I2_S12_0_MH2_MU2_0);
  t22 = t15 * t2;
  t23 = 0.2e1 * t22;
  t25 = 0.2e1 * t17 * t2;
  t26 = t1 * t1;
  t27 = RE(t26);
  t29 = 0.2e1 * t3 * t27;
  t30 = IM(t26);
  t32 = 0.2e1 * t9 * t30;
  t33 = IM(I2_MT2_MT2_MH2_MU2_0);
  t34 = t33 * t8;
  t35 = 0.2e1 * t34;
  t36 = t15 * s;
  t39 = 0.2e1 * t17 * s;
  t41 = t19 * t2 * s;
  t42 = t19 * t12;
  t43 = s * t22;
  t44 = IM(I2_S12_0_MH2_MU2_0);
  t46 = t44 * t8 * s;
  t47 = t34 * s;
  t48 = -0.8e1 * t19 * s - 0.16e2 * t3 * s - 0.16e2 * t10 + t14 - 0.16e2 * t15 + 0.16e2 * t17 + t23 - t25 - t29 + t32 - t35 + 0.6e1 * t36 - t39 + 0.16e2 * t4 + t41 + t42 - t43 - t46 + t47;
  t50 = beta_y * beta_y;
  t52 = t50 / s;
  t58 = 0.1e1 / (0.4e1 - s) / (beta_y - 0.1e1);
  t60 = 0.1e1 / (beta_y + 0.1e1);
  t65 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_2);
  t73 = RE(I3_MT2_MT2_S12_0_MT2_MH2_MU2_1);
  t77 = IM(I3_MT2_MT2_S12_0_MT2_MH2_MU2_1);
  t81 = 0.4e1 * s * t2 * t65 + 0.4e1 * s * t2 * t73 - 0.4e1 * s * t77 * t8 - 0.4e1 * t10 * s + 0.4e1 * s * t4 - t12 - t14 + t23 + t41 - t43 - t46 + t47;
  t94 = -0.2e1 * t12 * t65 - 0.2e1 * t12 * t73 - 0.2e1 * t27 * t65 - 0.2e1 * t27 * t73 + 0.2e1 * t30 * t77 + 0.4e1 * s - t25 - t29 + t32 - t35 + 0.2e1 * t36 + t39 - t42;
  return(0.16e2 * FH0 * t48 * t52 * t58 * t60 * PREF_V_CA_At - 0.32e2 * FA0 * (t81 + t94) * t52 * t58 * t60 * PREF_V_CA_Bt);
}
