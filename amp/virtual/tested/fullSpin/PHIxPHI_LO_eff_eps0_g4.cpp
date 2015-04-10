double Eval_B_EFF_PHIxPHI (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t10;
  double t11;
  double t14;
  double t15;
  double t16;
  double t2;
  double t24;
  double t26;
  double t38;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t1 = s * s;
  t2 = FA0 * FA0;
  t4 = FH0 * FH0;
  t5 = 0.4e1 * t2 + t4;
  t6 = t1 * t5;
  t7 = At * At;
  t8 = t7 * s;
  t9 = Bt * Bt;
  t10 = t9 * s;
  t11 = 0.4e1 * t7;
  t14 = DenS2(s, mH, GammaH);
  t15 = PREF_B_PHIxPHI * t14;
  t16 = sp(s1, s2);
  t24 = sp(s1, k2);
  t26 = sp(k1, s2);
  t38 = EPS_(k1, k2, s1, s2);
  return(0.8e1 * t6 * (-t8 + t10 + t11) * t15 * t16 + 0.16e2 * t6 * (At - Bt) * (At + Bt) * t24 * PREF_B_PHIxPHI * t14 * t26 - 0.8e1 * t6 * (-t8 - t10 + t11) * PREF_B_PHIxPHI * t14 + 0.32e2 * At * Bt * t1 * t5 * t38 * t15);
}


// double Eval_B_EFF_PHIxPHI (PS_2_2 const& ps)
// {
  
// #include "EVAL_V_PS_REFS.cpp"

//   std::cout << std::endl << " hallo! " << std::endl;
  
//   double t1;
//   double t10;
//   double t13;
//   double t14;
//   double t16;
//   double t17;
//   double t18;
//   double t19;
//   double t24;
//   double t28;
//   double t30;
//   double t33;
//   double t7;
//   double t8;
//   t1 = CA * CA;
//   t7 = VF(0.8e1 * CF * t1 * AlphaS2 / Pi2);
//   t8 = FA0 * FA0;
//   t10 = FH0 * FH0;
//   t13 = s * s;
//   t14 = DenS2(s, mH, GammaH);
//   t16 = At * At;
//   t17 = 1.0-4.0/s;
//   t18 = t16 * t17;
//   t19 = Bt * Bt;
//   t24 = sp(s2, s1);
//   t28 = sp(s2, k1);
//   t30 = sp(s1, k2);
//   t33 = EPS_(k1, k2, s1, s2);
//   return(t7 * (0.4e1 * t8 + t10) * t13 * t14 * (s * (t18 + t19) + s * (-t18 + t19) * t24 + (0.2e1 * t16 - 0.2e1 * t19) * t28 * t30 + 0.4e1 * At * Bt * t33));
// }
