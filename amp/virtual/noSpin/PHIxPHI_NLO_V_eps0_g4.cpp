
#include "AMP_HEADER.h"

#include <math.h>

double Eval_V_PHIxPHI (AMP_ARGS)
{
  using namespace Constants;
  AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_V(ap);

  double t1;
  double t11;
  double t12;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t23;
  double t25;
  double t3;
  double t31;
  double t33;
  double t39;
  double t4;
  double t54;
  double t56;
  double t7;
  double t8;
  t1 = s * s;
  t2 = t1 * s;
  t3 = beta * beta;
  t4 = t2 * t3;
  t7 = (RE(I2_S12_0_0_MU2_0)-2.0); // = log(RunParameters::MUR2 / s);
  t8 = t7 * t7;
  t11 = -CA * Pi2 + CA * t8 - 0.11e2;
  t12 = 0.1e1 / CA;
  t17 = RE(I1_MT2_MU2_0);
  t18 = t4 * t17;
  t20 = t1 * t1;
  t22 = t20 * (t3 + 0.1e1);
  t23 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t25 = t22 * t3 * t23;
  t31 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t33 = t20 * (beta - 0.1e1) * (beta + 0.1e1) * t3 * t31;
  t39 = t8 - Pi2 - 0.4e1;
  t54 = t2 * t17;
  t56 = t22 * t23;
  return(At2_fH2_De * (-0.4e1 * t4 * t11 * t12 * PREF_V_PHI_CA + (-0.8e1 * t18 - 0.4e1 * t25 + 0.8e1 * t33) * PREF_V_PHI_CF) + At2_fA2_De * (-0.16e2 * t4 * t39 * PREF_V_PHI_CA + (-0.32e2 * t18 - 0.16e2 * t25 + 0.32e2 * t33) * PREF_V_PHI_CF) + Bt2_fH2_De * (-0.4e1 * t2 * t11 * t12 * PREF_V_PHI_CA + (-0.8e1 * t54 - 0.4e1 * t56) * PREF_V_PHI_CF) + Bt2_fA2_De * (-0.16e2 * t2 * t39 * PREF_V_PHI_CA + (-0.32e2 * t54 - 0.16e2 * t56) * PREF_V_PHI_CF));
}
