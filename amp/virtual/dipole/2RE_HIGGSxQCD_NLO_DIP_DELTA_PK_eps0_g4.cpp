#include "AMP_HEADER.h"

double Eval_DIP_DELTA_PK (
 const double& s12_1,
 const double& t11_1,
 const double& t12_1,
 const double& B_1,
 const double& B_int_1,
 const double& Bt_1,
 const double& mu_F2,
 const double& mu_R2,
 const double& AlphaS)
{
  AMP_DEFINITIONS
    
  double t16;
  double t17;
  double t18;
  double t19;
  double t21;
  double t22;
  double t23;
  double t24;
  double t32;
  double t34;
  double t36;
  double t38;
  double t40;
  double t42;
  double t45;
  double t48;
  double t5;
  t5 = log(mu_F2 / s12_1);
  t16 = 0.1e1 / t11_1;
  t17 = 0.1e1 + t11_1;
  t18 = 0.1e1 / t17;
  t19 = log(t18);
  t21 = 0.1e1 / t12_1;
  t22 = 0.1e1 + t12_1;
  t23 = 0.1e1 / t22;
  t24 = log(t23);
  t32 = sqrt(t17);
  t34 = pow(0.1e1 - t32, 0.2e1);
  t36 = log(t34 * t16);
  t38 = sqrt(t22);
  t40 = pow(0.1e1 - t38, 0.2e1);
  t42 = log(t40 * t21);
  t45 = log(mu_F2 * t16);
  t48 = log(mu_F2 * t21);
  return(((-0.2e1 / 0.3e1 * Pi2 * CA - 0.4e1 * beta0 * t5) * B_1 + ((-0.67e2 / 0.18e2 + Pi2) * CA - 0.2e1 * beta0 + 0.5e1 / 0.9e1 * Nf) * B_int_1 + ((-t16 * t19 + t21 * t24 + (t11_1 - t12_1) * t23 * t18 / 0.2e1) * CA + (-0.2e1 * t36 + 0.2e1 * t42 + 0.2e1 * t45 - 0.2e1 * t48 + 0.4e1 * (-t38 + t32) / (t32 + 0.1e1) / (t38 + 0.1e1)) * beta0) * Bt_1) / TwoPi * AlphaS);
}
