double Eval_R_PHIxPHI_FSR (PS_2_3 const& ps)
{

#include "EVAL_R_PS_REFS.cpp"

  double t1;
  double t10;
  double t13;
  double t14;
  double t17;
  double t19;
  double t2;
  double t21;
  double t22;
  double t25;
  double t30;
  double t31;
  double t32;
  double t41;
  double t43;
  double t51;
  double t8;
  double t9;
  t1 = CF * CF;
  t2 = CA * CA;
  t8 = VF(0.128e3 * t1 * t2 * AlphaS3 / 0.3141592653589793e1);
  t9 = sp(p3, k1);
  t10 = t9 * t9;
  t13 = sp(p3, k2);
  t14 = t13 * t13;
  t17 = FA0 * FA0;
  t19 = FH0 * FH0;
  t21 = sp(p1, p2);
  t22 = t21 * t21;
  t25 = DenS2(p1 + p2, mH, GammaH);
  t30 = At * At;
  t31 = Bt * Bt;
  t32 = t30 + t31;
  t41 = t21 * t32;
  t43 = 0.4e1 * t30 - 0.2e1 * t41;
  t51 = 0.8e1 * t30;
  return(t8 / t10 / t14 * (0.4e1 * t17 + t19) * t22 * t25 * (t9 * t13 * (0.2e1 * t13 * ((0.4e1 * t13 - 0.8e1 * t21) * t32 + 0.16e2 * t30) + 0.4e1 * (0.2e1 - 0.2e1 * t21) * t43) + 0.2e1 * t10 * (0.2e1 * t13 * ((0.4e1 * t13 - 0.4e1 * t21) * t32 + t51) - 0.4e1 * t41 + t51) + 0.8e1 * t13 * t10 * t9 * t32 + 0.4e1 * t14 * t43) / 0.4e1);
}

