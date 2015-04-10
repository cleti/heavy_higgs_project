
#include "AMP_HEADER.h"

double Eval_V_2RE_PHI1xQCD0 (AMP_ARGS)
{
  using namespace Constants;
  AMP_DEFINITIONS
  HP_REFS_PHIxQCD(hp)
  AP_REFS_V(ap)

  double t1;
  double t10;
  double t15;
  double t16;
  double t17;
  double t18;
  double t2;
  double t23;
  double t24;
  double t26;
  double t27;
  double t29;
  double t35;
  double t37;
  double t38;
  double t39;
  double t4;
  double t45;
  double t55;
  double t59;
  double t6;
  double t68;
  double t7;
  double t73;
  double t75;
  double t78;
  double t9;
  t1 = CF * s;
  t2 = beta * beta;
  t4 = 0.1e1 / 0.3141592653589793e1;
  t6 = 0.1e1 / (beta_y - 0.1e1);
  t7 = t4 * t6;
  t9 = 0.1e1 / (beta_y + 0.1e1);
  t10 = RE(I1_MT2_MU2_0);
  t15 = s * s;
  t16 = CA * t15;
  t17 = t16 * t2;
  t18 = RE(I3_0_0_S12_0_0_0_MU2_0);
  t23 = CF * t15;
  t24 = t2 + 0.1e1;
  t26 = t23 * t2 * t24;
  t27 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t29 = t7 * t9 * t27;
  t35 = t15 * t2 * CF * (beta - 0.1e1);
  t37 = (beta + 0.1e1) * t4;
  t38 = t6 * t9;
  t39 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t45 = t7 * t9;
  t55 = t16 * t4;
  t59 = t23 * t24;
  t68 = IM(I3_0_0_S12_0_0_0_MU2_0);
  t73 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t75 = t7 * t9 * t73;
  t78 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  return(At_fH_re * (-0.32e2 * t1 * t10 * t2 * t7 * t9 - 0.32e2 * t17 * t18 * t7 * t9 + 0.32e2 * t35 * t37 * t38 * t39 + 0.176e3 * s * t2 * t45 - 0.16e2 * t26 * t29) * PREF_V + Bt_fA_re * (0.64e2 * t1 * t10 * t38 * t4 - 0.128e3 * CA * s * t45 + 0.64e2 * t18 * t38 * t55 + 0.32e2 * t29 * t59) * PREF_V + At_fH_im * (0.32e2 * t17 * t68 * t7 * t9 - 0.32e2 * t35 * t37 * t38 * t78 + 0.16e2 * t26 * t75) * PREF_V + Bt_fA_im * (-0.64e2 * t38 * t55 * t68 - 0.32e2 * t59 * t75) * PREF_V);
}
