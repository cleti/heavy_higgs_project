double Eval_V_PHIxPHI (
		       double const& s, 
		       double const& beta
		       )
{
  double t1;
  double t11;
  double t14;
  double t16;
  double t18;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t25;
  double t26;
  double t29;
  double t33;
  double t37;
  double t8;
  double t9;
  t1 = CF * CF;
  t2 = CA * CA;
  t8 = VF(0.8e1 * t1 * t2 * AlphaS3 / Pi3);
  t9 = FA0 * FA0;
  t11 = FH0 * FH0;
  t14 = s * s;
  t16 = DenS2(s, mH, GammaH);
  t18 = At * At;
  t19 = beta * beta;
  t20 = t19 * t19;
  t21 = t18 * t20;
  t22 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t23 = t22 * s;
  t25 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t26 = t25 * s;
  t29 = t18 * t19;
  t33 = Bt * Bt;
  t37 = RE(I1_MT2_MU2_0);
  return(-t8 * (0.4e1 * t9 + t11) * t14 * s * t16 * (s * t19 * t22 * t33 + s * t22 * t33 + t21 * t23 - 0.2e1 * t21 * t26 + t23 * t29 + 0.2e1 * t26 * t29 + 0.2e1 * t29 * t37 + 0.2e1 * t33 * t37) / 0.2e1);

}

