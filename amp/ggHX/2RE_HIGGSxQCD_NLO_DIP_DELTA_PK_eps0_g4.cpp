#include "AMP_HEADER.h"

double Eval_DIP_DELTA_PK (
  double const& s12_1,
  double const& t11_1,
  double const& t12_1,
  double const& B_1,
  double const& B_int_1,
  double const& Bt_1,
  double const& mu_F2,
  double const& mu_R2,
  double const& AlphaS)
{
  AMP_DEFINITIONS;
    
  double t5 = log(mu_F2 / s12_1);
  return((-0.2e1 / 0.3e1 * Pi2 * CA - 0.4e1 * beta0 * t5) * B_1 / TwoPi * AlphaS);
}
