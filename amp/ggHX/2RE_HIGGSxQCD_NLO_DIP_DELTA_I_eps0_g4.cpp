#include "AMP_HEADER.h"

double Eval_DIP_DELTA_I (
  double const& s12_1,
  double const& S12_1,
  double const& beta_1,
  double const& t11_1,
  double const& t12_1,
  double const& B_1,
  double const& Bt_1,
  double const& mu_R2)
{
  AMP_DEFINITIONS

  double t3 = log(MUR2 / s12_1);
  double t4 = 0.0;//t3 * t3; // cancels exactly against term in dV1
  return  ((0.67e2 / 0.9e1 - Pi2 + t4) * CA + 0.4e1 * beta0 * t3 + 0.4e1 * beta0 - 0.10e2 / 0.9e1 * Nf) * AlphaS / TwoPi * B_1;
}
