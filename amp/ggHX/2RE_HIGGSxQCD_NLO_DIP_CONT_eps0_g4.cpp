#include "AMP_HEADER.h"

double Eval_DIP_CONT (
  double const& s12_1,
  double const& t11_x,
  double const& t12_x,
  double const& B_x,
  double const& B_int_x,
  double const& Bt_x,
  double const& x,
  double const& mu_F2,
  double const& AlphaS)
{
  AMP_DEFINITIONS;
      
  // double t1 = 0.1e1 - x;
  // double t2 = 0.1e1 / t1;
  // double t7 = std::log(mu_F2 / (x*s12_x) * t2);
  // double t16 = P_gg_reg(x);
  // return (-((0.4e1 * t2 * t7) * CA + (0.2e1 * t7) * t16) * B_x * AlphaS / TwoPi);
double t2 = 0.1e1 / (0.1e1 - x);
double t8 = log(mu_F2 * t2 / x / s12_1);
double t12 = P_gg_reg(x);
return -(0.4e1 * t2 * t8 * CA + 0.2e1 * t8 * t12) * B_x * AlphaS / TwoPi;
}
