double Eval_R_PHIxPHI_FSR (PS_2_3 const& ps)
{

#include "EVAL_R_PS_REFS.cpp"

  double t1;
  double t10;
  double t100;
  double t101;
  double t104;
  double t110;
  double t119;
  double t124;
  double t128;
  double t13;
  double t14;
  double t17;
  double t19;
  double t2;
  double t21;
  double t22;
  double t25;
  double t26;
  double t30;
  double t31;
  double t32;
  double t41;
  double t43;
  double t51;
  double t67;
  double t68;
  double t72;
  double t74;
  double t75;
  double t8;
  double t9;
  double t91;
  double t93;
  double t94;
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
  t26 = t9 * t13;
  t30 = At * At;
  t31 = Bt * Bt;
  t32 = t30 + t31;
  t41 = t21 * t32;
  t43 = 0.4e1 * t30 - 0.2e1 * t41;
  t51 = 0.8e1 * t30;
  t67 = 0.2e1 * t9 + 0.2e1 * t13;
  t68 = t67 * t67;
  t72 = -0.8e1 * t13 * t21 * t9 + t68;
  t74 = sp(s2, k1);
  t75 = sp(s1, k2);
  t91 = sp(s2, s1);
  t93 = At * Bt;
  t94 = EPS_(k1, k2, s1, s2);
  t100 = t30 * (0.2e1 * t9 - 0.2e1 * t13);
  t101 = t31 * t67;
  t104 = sp(p3, s1);
  t110 = sp(p3, s2);
  t119 = EPS_(k1, p3, s1, s2);
  t124 = EPS_(k2, p3, s1, s2);
  t128 = t26 * (0.2e1 * t13 * ((0.4e1 * t13 - 0.8e1 * t21) * t32 + 0.16e2 * t30) + 0.4e1 * (0.2e1 - 0.2e1 * t21) * t43) + 0.2e1 * t10 * (0.2e1 * t13 * ((0.4e1 * t13 - 0.4e1 * t21) * t32 + t51) - 0.4e1 * t41 + t51) + 0.8e1 * t13 * t10 * t9 * t32 + 0.4e1 * t14 * t43 - 0.2e1 * (t30 - t31) * t72 * t74 * t75 - (t30 * (0.4e1 - 0.2e1 * t21) + 0.2e1 * t31 * t21) * (t68 + 0.4e1 * t26 * (0.2e1 * t9 + 0.2e1 * t13 - 0.2e1 * t21)) * t91 - 0.4e1 * t93 * t72 * t94 + 0.2e1 * t67 * (-t100 + t101) * t74 * t104 + 0.2e1 * t67 * (t100 + t101) * t75 * t110 + 0.2e1 * t32 * t68 * t104 * t110 - 0.4e1 * t93 * t9 * t67 * t119 + 0.4e1 * t93 * t13 * t67 * t124;
  return(t8 / t10 / t14 * (0.4e1 * t17 + t19) * t22 * t25 * t128 / 0.4e1);
}