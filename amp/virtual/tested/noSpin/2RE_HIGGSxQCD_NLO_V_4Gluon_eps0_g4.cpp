
#include "../../../../inc/Functions_Shared.h"

double Eval_V_4G (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  double t1;
  double t10;
  double t11;
  double t17;
  double t18;
  double t2;
  double t22;
  c_double t27;
  double t28;
  double t30;
  double t34;
  double t39;
  double t4;
  double t48;
  double t8;
  t1 = beta_y * beta_y;
  t2 = t1 * s;
  t4 = s * s;
  t8 = s * (0.4e1 * t2 - 0.3e1 * t4 + 0.16e2 * s - 0.16e2);
  t10 = 0.1e1 / (0.4e1 - s);
  t11 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t17 = s * (t2 - 0.2e1 * s + 0.8e1);
  t18 = RE(I2_MT2_0_MT2_MU2_0);
  t22 = RE(I2_S12_0_0_MU2_0);
  t27 = DenS(s, mH, GammaH);
  t28 = RE(t27);
  t30 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t34 = IM(I2_S12_0_0_MU2_0);
  t39 = IM(t27);
  t48 = VF(FH0 * CF * CA2 * AlphaS3 / 0.3141592653589793e1 * At);
  return(((-0.8e1 * t10 * t11 * t8 + 0.16e2 * t10 * t17 * t18 - 0.16e2 * t10 * t17 * t22) * t28 + (-0.16e2 * t10 * t17 * t34 - 0.8e1 * t10 * t30 * t8) * t39) * t48);
}
