
#include "AMP_HEADER.h"

double Eval_V_PHIxPHI (AMP_ARGS)
{

  using namespace Constants;
  AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_V(ap);

  double t1;
  double t10;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t27;
  double t28;
  double t30;
  double t31;
  double t34;
  double t39;
  double t4;
  double t40;
  double t41;
  double t42;
  double t45;
  double t46;
  double t49;
  double t5;
  double t51;
  double t52;
  double t6;
  double t9;
  t1 = s * s;
  t2 = t1 * s;
  t4 = beta * beta;
  t5 = CA * PREF_V_PHI_CF;
  t6 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t9 = t5 * t6 * t4 * s;
  t10 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t17 = (RE(I2_S12_0_0_MU2_0)-2.0); // =log(MUR2 / s);
  t18 = t17 * t17;
  t20 = CA * t18 * PREF_V_PHI_CA;
  t22 = t5 * t6 * s;
  t27 = CA * PREF_V_PHI_CA * Pi2;
  t28 = RE(I1_MT2_MU2_0);
  t30 = 0.2e1 * t5 * t28;
  t31 = 0.11e2 * PREF_V_PHI_CA;
  t34 = 0.1e1 / CA;
  t39 = PREF_V_PHI_CF * t6;
  t40 = t4 * s;
  t41 = t39 * t40;
  t42 = PREF_V_PHI_CF * t10;
  t45 = t18 * PREF_V_PHI_CA;
  t46 = t39 * s;
  t49 = PREF_V_PHI_CA * Pi2;
  t51 = 0.2e1 * PREF_V_PHI_CF * t28;
  t52 = 0.4e1 * PREF_V_PHI_CA;
  return(-0.4e1 * At2_fH2_De * t2 * t4 * (-0.2e1 * s * t10 * t4 * t5 + 0.2e1 * s * t10 * t5 + t20 + t22 - t27 + t30 - t31 + t9) * t34 - 0.16e2 * At2_fA2_De * t2 * t4 * (0.2e1 * s * t42 - 0.2e1 * t40 * t42 + t41 + t45 + t46 - t49 + t51 - t52) - 0.4e1 * Bt2_fH2_De * t2 * (t9 + t20 + t22 - t27 + t30 - t31) * t34 - 0.16e2 * Bt2_fA2_De * t2 * (t41 + t45 + t46 - t49 + t51 - t52));
}
