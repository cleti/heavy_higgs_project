
#include "AMP_HEADER.h"

double Eval_V_4G (AMP_ARGS)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_V(ap);

    

  double t1;
  double t10;
  double t14;
  double t19;
  double t23;
  double t24;
  double t3;
  double t31;
  double t35;
  double t4;
  double t7;
  double t9;
  t1 = beta * beta;
  t3 = y * y;
  t4 = t1 * s * t3;
  t7 = s * (t4 - 0.2e1 * s + 0.8e1);
  t9 = 0.1e1 / (0.4e1 - s);
  t10 = RE(I2_MT2_0_MT2_MU2_0);
  t14 = RE(I2_S12_0_0_MU2_0);
  t19 = s * s;
  t23 = s * (0.4e1 * t4 - 0.3e1 * t19 + 0.16e2 * s - 0.16e2);
  t24 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t31 = IM(I2_S12_0_0_MU2_0);
  t35 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  return(At_fH_re * (0.16e2 * t10 * t7 * t9 - 0.16e2 * t14 * t7 * t9 - 0.8e1 * t23 * t24 * t9) * PREF_V_CA + At_fH_im * (-0.8e1 * t23 * t35 * t9 - 0.16e2 * t31 * t7 * t9) * PREF_V_CA);
}
