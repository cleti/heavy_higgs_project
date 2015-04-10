
#include "../../../../inc/Functions_Shared.h"

double Eval_V_2RE_PHI1xQCD0 (
			     PS_2_2 const& ps
			     )
{
  #include "EVAL_V_PS_REFS"
  
  c_double t15;
  double t16;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t26;
  double t27;
  double t34;
  double t35;
  double t37;
  double t42;
  double t44;
  double t45;
  double t47;
  double t6;
  double t7;
  double t8;
  t6 = VF(0.16e2 * CF * CA * AlphaS3 / 0.3141592653589793e1);
  t7 = 1.0-4.0/s;
  t8 = beta_y * beta_y;
  t15 = DenS(s, mH, GammaH);
  t16 = IM(t15);
  t18 = IM(I3_0_0_S12_0_0_0_MU2_0);
  t19 = CA * s * t18;
  t20 = t7 + 0.1e1;
  t21 = t20 * s;
  t22 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t26 = (-t7 + 0.1e1) * s;
  t27 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t34 = RE(t15);
  t35 = RE(I3_0_0_S12_0_0_0_MU2_0);
  t37 = 0.2e1 * s * t35;
  t42 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t44 = t21 * t42 / 0.2e1;
  t45 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t47 = RE(I1_MT2_MU2_0);
  return(-t6 / (- t8 + 0.1e1) * (At * FH0 * s * t7 * (-t16 * (-0.2e1 * t19 + 0.2e1 * CF * (-t22 * t21 / 0.2e1 - t26 * t27)) + t34 * (-CA * (t37 - 0.11e2 / CA) + 0.2e1 * CF * (-t26 * t45 - t44 - t47))) - 0.2e1 * Bt * FA0 * s * (-t16 * (-CF * s * t20 * t22 - 0.2e1 * t19) + t34 * (-CA * (t37 - 0.4e1) + 0.2e1 * CF * (-t44 - t47)))));
}
