double Eval_B_QCDxQCD_GG (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"
  
  double t1;
  double t11;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t28;
  double t3;
  double t30;
  double t31;
  double t37;
  double t38;
  double t4;
  double t46;
  double t48;
  double t5;
  double t53;
  double t54;
  double t55;
  double t56;
  double t6;
  double t7;
  double t8;
  t1 = beta_y * beta_y;
  t2 = t1 * t1;
  t3 = s * s;
  t4 = t2 * t3;
  t5 = sp(s2, p1);
  t6 = sp(s1, p1);
  t7 = t5 * t6;
  t8 = t1 * s;
  t11 = 0.8e1 * t8;
  t16 = 0.8e1 * s;
  t18 = 0.1e1 / t3;
  t19 = (-0.8e1 * s * t5 * t6 + 0.8e1 * t7 * t8 + t11 - t16 - t3 + t4 + 0.32e2 * t7 + 0.32e2) * t18;
  t20 = beta_y + 0.1e1;
  t21 = 0.1e1 / t20;
  t22 = beta_y - 0.1e1;
  t23 = 0.1e1 / t22;
  t28 = t8 - s + 0.4e1;
  t30 = t28 * t5 * t18;
  t31 = sp(s1, k2);
  t37 = t28 * t6 * t18;
  t38 = sp(k1, s2);
  t46 = (-0.2e1 * t1 * t3 + t11 - t16 + t3 + t4 + 0.32e2) * t18;
  t48 = sp(s1, s2);
  t53 = t20 * t20;
  t54 = 0.1e1 / t53;
  t55 = t22 * t22;
  t56 = 0.1e1 / t55;
  return(0.16e2 * t21 * t23 * t46 * t48 * PREF_B_QCDxQCD_CA + 0.256e3 * t21 * t37 * t38 * t56 * PREF_B_QCDxQCD_CF - 0.256e3 * t23 * t30 * t31 * t54 * PREF_B_QCDxQCD_CF + 0.64e2 * t46 * t48 * t54 * t56 * PREF_B_QCDxQCD_CF - 0.16e2 * t19 * t21 * t23 * PREF_B_QCDxQCD_CA - 0.64e2 * t19 * t54 * t56 * PREF_B_QCDxQCD_CF - 0.64e2 * t21 * t30 * t31 * PREF_B_QCDxQCD_CA + 0.64e2 * t23 * t37 * t38 * PREF_B_QCDxQCD_CA);
}
