double Eval_R_PHIxPHI_ISR (PS_2_3 const& ps)
{

#include "EVAL_R_PS_REFS.cpp"

    double t1;
  double t11;
  double t14;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t23;
  double t24;
  double t25;
  double t27;
  double t28;
  double t29;
  double t39;
  double t40;
  double t41;
  double t42;
  double t47;
  double t51;
  double t53;
  double t56;
  double t8;
  double t9;
  t1 = CA * CA;
  t8 = VF(0.64e2 * CF * t1 * CA * AlphaS3 / 0.3141592653589793e1);
  t9 = FA0 * FA0;
  t11 = FH0 * FH0;
  t14 = sp(k1, k2);
  t17 = pow(0.2e1 + 0.2e1 * t14, 0.2e1);
  t18 = t17 * t17;
  t19 = sp(p1, p2);
  t20 = t19 * t19;
  t21 = t20 * t20;
  t23 = sp(p1, p3);
  t24 = t23 * t23;
  t25 = t24 * t24;
  t27 = sp(p3, p2);
  t28 = t27 * t27;
  t29 = t28 * t28;
  t39 = DenS2(k1 + k2, mH, GammaH);
  t40 = At * At;
  t41 = t40 * beta_tt2;
  t42 = Bt * Bt;
  t47 = sp(s2, s1);
  t51 = sp(s2, k1);
  t53 = sp(s1, k2);
  t56 = EPS_(k1, k2, s1, s2);
  return(t8 * (0.4e1 * t9 + t11) * (t18 + 0.16e2 * t21 + 0.16e2 * t25 + 0.16e2 * t29) / t19 / t23 / t27 * t39 * (s_tt * (t41 + t42) + s_tt * (-t41 + t42) * t47 + (0.2e1 * t40 - 0.2e1 * t42) * t51 * t53 + 0.4e1 * At * Bt * t56) / 0.8e1);
}
