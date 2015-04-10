#include <math.h>

double Eval_DIP_CONT (
  double s12_x,
  double t11_x,
  double t12_x,
  double B_x,
  double B_int_x,
  double Bt_x,
  double x,
  double mu_F2)
{
  double t1;
  double t16;
  double t2;
  double t21;
  double t28;
  double t3;
  double t30;
  double t34;
  double t36;
  double t39;
  double t43;
  double t45;
  double t48;
  double t53;
  double t58;
  double t7;
  double t9;
  t1 = 0.1e1 - x;
  t2 = 0.1e1 / t1;
  t3 = 0.1e1 / x;
  t7 = log(mu_F2 * t3 / s12_x);
  t9 = log(t1);
  t16 = P_gg_reg(x);
  t21 = log(t1 * t3);
  t28 = 0.1e1 / t12_x;
  t30 = log(t11_x * t28);
  t34 = t1 * t11_x + 0.1e1;
  t36 = t1 * t12_x + 0.1e1;
  t39 = log(t34 / t36);
  t43 = log(t34);
  t45 = 0.1e1 / t11_x;
  t48 = log(t36);
  t53 = pow(0.1e1 - x + t28, 0.2e1);
  t58 = pow(0.1e1 - x + t45, 0.2e1);
  return(AlphaS * (((-0.4e1 * t2 * t7 + 0.4e1 * t2 * t9) * CA + (0.2e1 * t9 - 0.2e1 * t7) * t16) * B_x + (0.2e1 * CA * t2 * t21 + t16 * t21) * B_int_x + ((-0.2e1 * t30 * t2 + 0.2e1 * (t39 - t30) * t2 + 0.2e1 * t43 * t3 * t45 - 0.2e1 * t48 * t3 * t28 + t1 / t53 / 0.2e1 - t1 / t58 / 0.2e1) * CA + (t39 - 0.2e1 * t30) * t16) * Bt_x) / TwoPi);
}
