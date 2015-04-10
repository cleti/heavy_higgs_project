double Eval_UID_PHIxPHI_SE00  (
		      PS_2_3 const& ps,
		      PS_2_2 const& ps_red
		      ) 
{
#include "EVAL_UID_PS_REFS.cpp"

  double t1;
  double t10;
  double t11;
  double t12;
  double t13;
  double t16;
  double t17;
  double t18;
  double t2;
  double t20;
  double t22;
  double t26;
  double t27;
  double t28;
  double t3;
  double t30;
  double t34;
  double t36;
  double t38;
  double t4;
  double t43;
  double t47;
  double t5;
  double t7;
  double t8;
  t1 = sp(P2, p1);
  t2 = sp(P2, p3);
  t3 = t1 * t2;
  t4 = sp(p3, p1);
  t5 = t3 * t4;
  t7 = DenS2(P2 + p1, mH, GammaH);
  t8 = FH0 * FH0;
  t10 = CA * CA;
  t11 = t10 * CA;
  t12 = t8 * CF * t11;
  t13 = AlphaS * AlphaS;
  t16 = t13 * AlphaS / 0.3141592653589793e1;
  t17 = Bt * Bt;
  t18 = t17 * t16;
  t20 = VF(t12 * t18);
  t22 = sp(K2, K1);
  t26 = t7 * t22;
  t27 = At * At;
  t28 = t16 * t27;
  t30 = VF(t12 * t28);
  t34 = FA0 * FA0;
  t36 = t34 * CF * t11;
  t38 = VF(t36 * t18);
  t43 = VF(t36 * t28);
  t47 = t4 * t7;
  return(0.1024e4 * t20 * t22 * t5 * t7 + 0.1024e4 * t20 * t3 * t47 + 0.1024e4 * t26 * t30 * t5 + 0.4096e4 * t26 * t38 * t5 + 0.4096e4 * t26 * t43 * t5 - 0.1024e4 * t3 * t30 * t47 + 0.4096e4 * t3 * t38 * t47 - 0.4096e4 * t3 * t43 * t47);
}
