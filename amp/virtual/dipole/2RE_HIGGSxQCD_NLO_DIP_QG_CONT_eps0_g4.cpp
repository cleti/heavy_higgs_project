

double Eval_DIP_QG_CONT (
			 double const& s12_1,
			 double const& t11_x,
			 double const& t12_x,
			 double const& B_x,
			 double const& Bt_x,
			 double const& x,
			 double const& mu_F2,
			 double const& AlphaS)
{
  using namespace Constants;
  double t1;
  double t10;
  double t14;
  double t15;
  double t16;
  double t18;
  double t2;
  double t21;
  double t22;
  double t3;
  double t30;
  double t34;
  double t7;
  t1 = P_qg_reg(x);
  t2 = 0.1e1 - x;
  t3 = log(t2);
  t7 = log(mu_F2 / s12_1);
  t10 = P1_qg(x);
  t14 = t11_x * t2 + 0.1e1;
  t15 = log(t14);
  t16 = 0.1e1 / x;
  t18 = 0.1e1 / t11_x;
  t21 = t12_x * t2 + 0.1e1;
  t22 = log(t21);
  t30 = log(t12_x * t18);
  t34 = log(t14 / t21);
  return(AlphaS * ((t1 * (0.2e1 * t3 - t7) + t10) * B_x + (0.2e1 * (t15 * t16 * t18 - t22 * t16 / t12_x) * CF + t1 * (0.2e1 * t30 + t34)) * Bt_x) / TwoPi);
}
