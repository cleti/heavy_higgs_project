double Eval_B_QCDxQCD_QQ (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t12;
  double t13;
  double t19;
  double t2;
  double t20;
  double t6;
  double t7;
  double t9;
  t1 = 0.1e1 - 0.4e1 / s;
  t2 = beta_y * beta_y;
  t6 = t1 * (-t2 / t1 + 0.1e1);
  t7 = sp(s1, s2);
  t9 = 0.1e1 / s;
  t12 = sp(s1, p1);
  t13 = sp(s2, p2);
  t19 = sp(s1, p2);
  t20 = sp(s2, p1);
  return(0.8e1 * PREF_B_QCDxQCD * (0.2e1 - t6 + t6 * t7 - 0.4e1 * t9 * (beta_y + 0.1e1) * t12 * t13 - 0.4e1 * t9 * (-beta_y + 0.1e1) * t19 * t20));
}
