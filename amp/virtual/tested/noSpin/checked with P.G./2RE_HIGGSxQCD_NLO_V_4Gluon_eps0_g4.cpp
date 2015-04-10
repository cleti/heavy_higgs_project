double Eval_V_4G (
		  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t10;
  double t12;
  double t13;
  double t19;
  double t20;
  double t24;
  c_double t29;
  double t30;
  double t32;
  double t36;
  double t4;
  double t41;
  double t50;
  double t6;
  t1 = beta_y * beta_y;
  t4 = t1 * s;
  t6 = s * s;
  t10 = s * (0.4e1 * t4 - 0.3e1 * t6 + 0.16e2 * s - 0.16e2);
  t12 = 0.1e1 / (0.4e1 - s);
  t13 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t19 = s * (t4 - 0.2e1 * s + 0.8e1);
  t20 = RE(I2_MT2_0_MT2_MU2_0);
  t24 = RE(I2_S12_0_0_MU2_0);
  t29 = DenS(s, mH, GammaH);
  t30 = RE(t29);
  t32 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t36 = IM(I2_S12_0_0_MU2_0);
  t41 = IM(t29);
  t50 = VF(FH0 * CF * CA2 * AlphaS3 / 0.3141592653589793e1 * At);
  return(((-0.8e1 * t10 * t12 * t13 + 0.16e2 * t12 * t19 * t20 - 0.16e2 * t12 * t19 * t24) * t30 + (-0.8e1 * t10 * t12 * t32 - 0.16e2 * t12 * t19 * t36) * t41) * t50);
}
