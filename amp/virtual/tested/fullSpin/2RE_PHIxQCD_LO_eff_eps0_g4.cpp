double Eval_B_EFF_2PHIxQCD (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t12;
  c_double t14;
  double t15;
  double t17;
  double t2;
  double t25;
  double t26;
  double t27;
  double t3;
  double t39;
  double t40;
  double t44;
  double t5;
  double t51;
  double t52;
  double t6;
  double t9;
  t1 = At * FH0;
  t2 = t1 * s;
  t3 = Bt * FA0;
  t5 = 0.2e1 * t3 * s;
  t6 = 0.4e1 * t1;
  t9 = 0.1e1 / (beta_y - 0.1e1);
  t12 = 0.1e1 / (beta_y + 0.1e1);
  t14 = DenS(s, mH, GammaH);
  t15 = RE(t14);
  t17 = sp(s1, s2);
  t25 = sp(s1, k2);
  t26 = t25 * PREF_B_PHIxQCD;
  t27 = sp(k1, s2);
  t39 = 0.2e1 * At * FA0;
  t40 = Bt * FH0;
  t44 = EPS_(k1, k2, s1, s2);
  t51 = (t39 + t40) * t9 * t12;
  t52 = IM(t14);
  return(0.64e2 * (-t2 - t5 + t6) * t9 * t12 * PREF_B_PHIxQCD * t15 * t17 + 0.128e3 * (t1 + 0.2e1 * t3) * t9 * t12 * t26 * t15 * t27 - 0.64e2 * (-t2 + t5 + t6) * t9 * t12 * PREF_B_PHIxQCD * t15 - 0.128e3 * (t39 - t40) * t9 * t12 * t44 * PREF_B_PHIxQCD * t15 - 0.128e3 * t51 * PREF_B_PHIxQCD * t52 * t27 - 0.128e3 * t51 * t26 * t52);
}
