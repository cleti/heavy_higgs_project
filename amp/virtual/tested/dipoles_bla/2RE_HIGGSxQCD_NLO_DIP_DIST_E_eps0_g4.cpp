#include <math.h>

double Eval_DIP_DIST_E (
  double s12_1,
  double t11_1,
  double t12_1,
  double B_1,
  double B_int_1,
  double Bt_1,
  double x,
  double mu_F2)
{
  double t1;
  double t15;
  double t2;
  double t20;
  double t22;
  double t31;
  double t37;
  double t42;
  double t5;
  double t7;
  t1 = 0.1e1 - x;
  t2 = 0.1e1 / t1;
  t5 = log(mu_F2 / s12_1);
  t7 = log(t1);
  t15 = log(t1 / x);
  t20 = 0.1e1 / t12_1;
  t22 = log(t11_1 * t20);
  t31 = log((t1 * t11_1 + 0.1e1) / (t1 * t12_1 + 0.1e1));
  t37 = pow(0.1e1 - x + 0.1e1 / t11_1, 0.2e1);
  t42 = pow(0.1e1 - x + t20, 0.2e1);
  return(-AlphaS * ((-0.4e1 * t2 * t5 + 0.4e1 * t2 * t7) * CA * B_1 + 0.2e1 * CA * t15 * t2 * B_int_1 + (-0.2e1 * t22 * t2 + 0.2e1 * (t31 - t22) * t2 - t1 / t37 / 0.2e1 + t1 / t42 / 0.2e1) * CA * Bt_1) / TwoPi);
}
