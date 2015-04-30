

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
  double t13;
  double t3;
  double t7;
  t3 = P_qg_reg(x);
  t7 = log(mu_F2 / s12_1 / x );
  t13 = log(t12_x / t11_x);
  return(AlphaS / TwoPi * (t3 * (- t7) * B_x + t3 * t13 * Bt_x));
}
