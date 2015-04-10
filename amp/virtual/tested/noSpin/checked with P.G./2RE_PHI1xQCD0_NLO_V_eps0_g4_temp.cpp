double Eval_V_2RE_PHI1xQCD0 (double s, double y, double beta)
{
  c_double t15;
  double t16;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t26;
  double t27;
  double t34;
  double t35;
  double t37;
  double t40;
  double t43;
  double t45;
  double t46;
  double t48;
  double t6;
  double t7;
  double t8;
  t6 = VF(0.16e2 * CF * CA * AlphaS3 / 0.3141592653589793e1);
  t7 = beta * beta;
  t8 = y * y;
  t15 = DenS(s, mH, GammaH);
  t16 = IM(t15);
  t18 = IM(I3_0_0_S12_0_0_0_MU2_0);
  t19 = CA * s * t18;
  t20 = t7 + 0.1e1;
  t21 = t20 * s;
  t22 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t26 = (-t7 + 0.1e1) * s;
  t27 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t34 = RE(t15);
  t35 = RE(I3_0_0_S12_0_0_0_MU2_0);
  t37 = 0.2e1 * s * t35;
  t40 = 0.11e2 / CA;
  t43 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t45 = t21 * t43 / 0.2e1;
  t46 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t48 = RE(I1_MT2_MU2_0);
  return(-t6 / (-t7 * t8 + 0.1e1) * (At * FH0 * s * t7 * (t16 * (-0.2e1 * t19 + 0.2e1 * CF * (-t22 * t21 / 0.2e1 - t26 * t27)) + t34 * (-CA * (t37 - t40) + 0.2e1 * CF * (-t26 * t46 - t45 - t48))) - 0.2e1 * Bt * FA0 * s * (t16 * (-CF * s * t20 * t22 - 0.2e1 * t19) + t34 * (-CA * (t37 - t40 - 0.4e1) + 0.2e1 * CF * (-t45 - t48)))));
}
