double Eval_UID_00SE (
		      PS_2_3 const& ps,
		      PS_2_2 const& ps_red
		      ) 
{
#include "EVAL_UID_PS_REFS.cpp"

  double t1;
  double t10;
  double t11;
  double t14;
  double t15;
  double t17;
  double t24;
  double t25;
  double t26;
  double t29;
  double t3;
  double t31;
  double t32;
  double t35;
  double t4;
  double t40;
  double t41;
  double t46;
  double t47;
  double t51;
  double t57;
  double t6;
  double t7;
  double t71;
  c_double t77;
  double t78;
  double t80;
  double t9;
  t1 = CF * CF;
  t3 = FA0 * t1 * CA;
  t4 = AlphaS * AlphaS;
  t6 = t4 * AlphaS * 0.3141592653589793e1;
  t7 = t6 * At;
  t9 = VF(t3 * t7);
  t10 = sp(K1, p1);
  t11 = sp(K1, p2);
  t14 = 0.1e1 / t10;
  t15 = 0.1e1 / t11;
  t17 = EPS_(K1, K2, p1, p2);
  t24 = VF(FH0 * t1 * CA * t7);
  t25 = t10 * t10;
  t26 = t24 * t25;
  t29 = sp(K2, p2);
  t31 = t24 * t10;
  t32 = t11 * t11;
  t35 = sp(K2, p1);
  t40 = sp(p2, p1);
  t41 = t11 * t40;
  t46 = sp(K2, K1);
  t47 = t40 * t46;
  t51 = t24 * t11;
  t57 = VF(t3 * t6 * Bt);
  t71 = -0.4e1 * t10 * t11 * t40 * t57 - 0.2e1 * t10 * t29 * t40 * t57 - t11 * t29 * t31 - t11 * t31 * t35 - t24 * t32 * t35 + t29 * t31 * t40 + t35 * t40 * t51 - 0.2e1 * t35 * t41 * t57 + 0.2e1 * t11 * t26 - t26 * t29 + 0.2e1 * t31 * t32 - t31 * t40 - 0.2e1 * t31 * t41 + t31 * t47 - t40 * t51 + t47 * t51;
  t77 = DenS(p1 + p2, mH, GammaH);
  t78 = RE(t77);
  t80 = VQgFF(k2, p3, k1);
  return((0.512e3 * t9 * (t10 - t11) * t14 * t15 * t17 - 0.256e3 * t71 * t14 * t15) * t78 * t80);
}
