

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
  double t12;
  double t13;
  double t14;
  double t16;
  double t19;
  double t2;
  double t20;
  double t28;
  double t3;
  double t31;
  double t5;
  double t8;
  t1 = P_qg_reg(x);
  t2 = 0.1e1 - x;
  t3 = log(t2);
  t5 = log(x);
  t8 = P1_qg(x);
  t12 = t11_x * t2 + 0.1e1;
  t13 = log(t12);
  t14 = 0.1e1 / x;
  t16 = 0.1e1 / t11_x;
  t19 = t12_x * t2 + 0.1e1;
  t20 = log(t19);
  t28 = log(t12_x * t16);
  t31 = log(t12 / t19);
  return(AlphaS * ((t1 * (0.2e1 * t3 - t5) + t8) * B_x + (0.2e1 * (t13 * t14 * t16 - t20 * t14 / t12_x) * CF + t1 * (t28 + t31)) * Bt_x) / TwoPi);
}

