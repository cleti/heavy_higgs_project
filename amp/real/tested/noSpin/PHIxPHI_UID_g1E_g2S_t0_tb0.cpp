double Eval_UID_PHIxPHI_ES00 (
		      PS_2_3 const& ps,
		      PS_2_2 const& ps_red
		      ) 
{
#include "EVAL_UID_PS_REFS.cpp"

  double t1;
  double t10;
  double t11;
  double t13;
  double t14;
  double t16;
  double t18;
  double t19;
  double t22;
  double t23;
  double t26;
  double t27;
  double t29;
  double t3;
  double t34;
  double t36;
  double t38;
  double t4;
  double t44;
  double t5;
  double t6;
  double t9;
  t1 = FH0 * FH0;
  t3 = CA * CA;
  t4 = t3 * CA;
  t5 = t1 * CF * t4;
  t6 = AlphaS * AlphaS;
  t9 = t6 * AlphaS / 0.3141592653589793e1;
  t10 = Bt * Bt;
  t11 = t9 * t10;
  t13 = VF(t5 * t11);
  t14 = sp(K2, K1);
  t16 = sp(P1, p2);
  t18 = sp(P1, p3);
  t19 = sp(p3, p2);
  t22 = DenS2(P1 + p2, mH, GammaH);
  t23 = t18 * t19 * t22;
  t26 = At * At;
  t27 = t9 * t26;
  t29 = VF(t5 * t27);
  t34 = FA0 * FA0;
  t36 = t34 * CF * t4;
  t38 = VF(t36 * t11);
  t44 = VF(t36 * t27);
  return(0.1024e4 * t13 * t14 * t16 * t23 + 0.1024e4 * t14 * t16 * t23 * t29 + 0.4096e4 * t14 * t16 * t23 * t38 + 0.4096e4 * t14 * t16 * t23 * t44 + 0.1024e4 * t13 * t16 * t23 - 0.1024e4 * t16 * t23 * t29 + 0.4096e4 * t16 * t23 * t38 - 0.4096e4 * t16 * t23 * t44);
}
