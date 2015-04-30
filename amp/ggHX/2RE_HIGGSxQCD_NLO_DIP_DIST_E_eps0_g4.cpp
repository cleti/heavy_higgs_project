#include "AMP_HEADER.h"

double Eval_DIP_DIST_E (
  double const& s12_1,
  double const& t11_1,
  double const& t12_1,
  double const& B_1,
  double const& B_int_1,
  double const& Bt_1,
  double const& x,
  double const& mu_F2,
  double const& AlphaS)
{
  AMP_DEFINITIONS;
    
  // double t1 = 0.1e1 - x;
  // double t2 = 0.1e1 / t1;
  // double t5 = std::log(mu_F2 / s12_1 * t2);
  // return((-0.4e1 * t2 * t5) * CA * B_1 * AlphaS / TwoPi);
  double t2 = 0.1e1 / (0.1e1 - x);
  double t6 = log(mu_F2 * t2 / s12_1);
  return -0.4e1 * t2 * t6 * CA * B_1 * AlphaS / TwoPi;
}
