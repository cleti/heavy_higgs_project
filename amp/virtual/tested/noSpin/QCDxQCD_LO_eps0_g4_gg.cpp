#include <math.h>

double Eval_B_QCDxQCD_GG (PS_2_2 const& ps)
{
  
#include "EVAL_V_PS_REFS"

  double t1;
  double t10;
  double t13;
  double t16;
  double t17;
  double t19;
  double t22;
  double t27;
  double t28;
  double t30;
  double t37;
  double t39;
  double t46;
  double t5;
  double t7;
  t1 = beta_y * beta_y;
  t5 = Pi2;
  t7 = VF(CF * AlphaS2 * t5);
  t10 = 1-4.0/s;
  t13 = -t1 / t10 + 0.1e1;
  t16 = t10 * t10;
  t17 = t13 * t13;
  t19 = t16 * (t17 + 0.1e1);
  t22 = 0.0;//sp(s1, s2);
  t27 = 0.1e1 / s;
  t28 = 0.0;//sp(s1, p1);
  t30 = 0.0;//sp(s2, p2);
  t37 = 0.0;//sp(s1, p2);
  t39 = 0.0;//sp(s2, p1);
  t46 = pow(-t1 + 0.1e1, 0.2e1);
  return(0.16e2 * (PREF_B_QCDxQCD_CA * (t1 + 0.1e1) - 0.2e1 * t7) * (0.1e1 + 0.2e1 * t10 * t13 - t19 + (0.1e1 - 0.2e1 * t10 + t19) * t22 - 0.4e1 * (beta_y + 0.1e1) * t10 * t13 * t27 * t28 * t30 - 0.4e1 * (-beta_y + 0.1e1) * t10 * t13 * t27 * t37 * t39) / t46);
}
