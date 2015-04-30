
#include "AMP_HEADER.h"

double Eval_V_PHIxPHI_IM_INTab (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_V(ap);
  
  // need the imaginary parts of these prefactors here !!!
  double const& At_Bt_fH2_De = hp.At_Bt_fH2_DeIM;
  double const& At_Bt_fA2_De = hp.At_Bt_fA2_DeIM;


  double t1;
  double t10;
  double t106;
  double t11;
  double t14;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t23;
  double t25;
  double t27;
  double t29;
  double t31;
  double t32;
  double t35;
  double t4;
  double t48;
  double t49;
  double t50;
  double t53;
  double t55;
  double t59;
  double t6;
  double t61;
  double t74;
  double t78;
  double t8;
  double t85;
  double t86;
  double t88;
  double t95;
  t1 = s * s;
  t2 = sp(k1, s2);
  t4 = t2 * At_Bt_fH2_De * PREF_V_PHI_CA;
  t6 = sp(s1, k2);
  t8 = t6 * At_Bt_fH2_De * PREF_V_PHI_CA;
  t10 = CA * t2;
  t11 = At_Bt_fA2_De * PREF_V_PHI_CA;
  t14 = CA * t6;
  t17 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t18 = t10 * t17;
  t20 = s * At_Bt_fH2_De * PREF_V_PHI_CF;
  t22 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t23 = t10 * t22;
  t25 = t14 * t17;
  t27 = t14 * t22;
  t29 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t31 = EPS_(k1, k2, s1, s2);
  t32 = CA * t29 * t31;
  t35 = s * At_Bt_fA2_De * PREF_V_PHI_CF;
  t48 = (RE(I2_S12_0_0_MU2_0)-2.0); // = log(RunParameters::MUR2 / s);
  t49 = t48 * t48;
  t50 = CA * t49;
  t53 = RE(I1_MT2_MU2_0);
  t55 = t53 * At_Bt_fH2_De * PREF_V_PHI_CF;
  t59 = Pi2 * At_Bt_fH2_De * PREF_V_PHI_CA;
  t61 = -0.16e2 * t10 * t11 + 0.2e1 * t10 * t55 - t10 * t59 - 0.16e2 * t11 * t14 + t18 * t20 + 0.4e1 * t18 * t35 + t20 * t23 + t20 * t25 + t20 * t27 - t20 * t32 + 0.4e1 * t23 * t35 + 0.4e1 * t25 * t35 + 0.4e1 * t27 * t35 - 0.4e1 * t32 * t35 + t4 * t50 + t50 * t8 - 0.11e2 * t4 - 0.11e2 * t8;
  t74 = t53 * At_Bt_fA2_De * PREF_V_PHI_CF;
  t78 = Pi2 * At_Bt_fA2_De * PREF_V_PHI_CA;
  t85 = beta * beta;
  t86 = t85 * s;
  t88 = t86 * At_Bt_fH2_De * PREF_V_PHI_CF;
  t95 = t86 * At_Bt_fA2_De * PREF_V_PHI_CF;
  t106 = 0.4e1 * t2 * t50 * At_Bt_fA2_De * PREF_V_PHI_CA + 0.4e1 * t50 * t6 * At_Bt_fA2_De * PREF_V_PHI_CA + 0.8e1 * t10 * t74 - 0.4e1 * t10 * t78 + 0.2e1 * t14 * t55 - t14 * t59 + 0.8e1 * t14 * t74 - 0.4e1 * t14 * t78 + t18 * t88 + 0.4e1 * t18 * t95 - t23 * t88 - 0.4e1 * t23 * t95 + t25 * t88 + 0.4e1 * t25 * t95 - t27 * t88 - 0.4e1 * t27 * t95 + t32 * t88 + 0.4e1 * t32 * t95;
  return(0.16e2 * t1 * (t61 + t106) / CA);
}
