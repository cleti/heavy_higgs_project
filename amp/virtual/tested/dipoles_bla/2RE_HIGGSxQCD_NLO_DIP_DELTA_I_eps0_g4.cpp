#include <math.h>

double Eval_DIP_DELTA_I (
  double s12_1,
  double S12_1,
  double beta_1,
  double t11_1,
  double t12_1,
  double B_1,
  double Bt_1,
  double mu_R2)
{
  double rho;
  double t100;
  double t101;
  double t104;
  double t105;
  double t106;
  double t107;
  double t108;
  double t112;
  double t115;
  double t117;
  double t120;
  double t125;
  double t127;
  double t129;
  double t13;
  double t135;
  double t139;
  double t141;
  double t145;
  double t18;
  double t2;
  double t21;
  double t25;
  double t3;
  double t32;
  double t34;
  double t35;
  double t38;
  double t39;
  double t4;
  double t42;
  double t46;
  double t5;
  double t50;
  double t53;
  double t54;
  double t59;
  double t6;
  double t7;
  double t73;
  double t77;
  double t78;
  double t8;
  double t87;
  double t88;
  double t89;
  double t9;
  double t92;
  double t93;
  double t94;
  double t95;
  double t96;
  double t99;
  t2 = 0.1e1 + beta_1;
  t3 = 0.1e1 / t2;
  rho = (0.1e1 - beta_1) * t3;
  t4 = beta_1 * beta_1;
  t5 = t4 + 0.1e1;
  t6 = 0.1e1 / beta_1;
  t7 = t5 * t6;
  t8 = log(rho);
  t9 = t8 * t8;
  t13 = log(t5 / 0.2e1);
  t18 = log(mu_R2 / S12_1);
  t21 = t2 * t2;
  t25 = log(0.2e1 * t5 / t21);
  t32 = log(0.1e1 / mu_R2);
  t34 = 0.1e1 / t5;
  t35 = t4 - 0.1e1;
  t38 = sqrt(-t35);
  t39 = 0.2e1 - t38;
  t42 = log(t38 / t39);
  t46 = log(0.1e1 - t38);
  t50 = log(0.1e1 - t38 / 0.2e1);
  t53 = rho * rho;
  t54 = dilog(t53);
  t59 = dilog(0.2e1 * beta_1 * t3);
  t73 = -t7 * t9 / 0.2e1 + (-t13 * t5 * t6 + t18 * t5 * t6 + 0.2e1 * t25 * t5 * t6) * t8 + t32 + 0.3e1 * t18 + 0.2e1 * t34 * t35 * t42 - 0.4e1 * t46 + 0.2e1 * t50 + 0.3e1 * t13 + 0.2e1 * t7 * t54 - 0.2e1 * t7 * t59 - t5 * Pi2 * t6 / 0.3e1 - 0.4e1 * (t38 * t4 + 0.5e1 * t38 - 0.6e1) / t39 * t34;
  t77 = log(mu_R2 / s12_1);
  t78 = t77 * t77;
  t87 = 0.1e1 / t11_1;
  t88 = log(t87);
  t89 = t88 * t88;
  t92 = log(mu_R2 * t87);
  t93 = 0.1e1 + t11_1;
  t94 = 0.1e1 / t93;
  t95 = t11_1 * t94;
  t96 = log(t95);
  t99 = 0.1e1 / t12_1;
  t100 = log(t99);
  t101 = t100 * t100;
  t104 = log(mu_R2 * t99);
  t105 = 0.1e1 + t12_1;
  t106 = 0.1e1 / t105;
  t107 = t12_1 * t106;
  t108 = log(t107);
  t112 = log(t94);
  t115 = log(t106);
  t117 = t92 * t92;
  t120 = t104 * t104;
  t125 = dilog(t95);
  t127 = dilog(t107);
  t129 = t89 / 0.2e1 + (-t92 + t96) * t88 - t101 / 0.2e1 + (t104 - t108) * t100 + (t96 + t87) * t112 + (-t108 - t99) * t115 - t117 / 0.2e1 - 0.3e1 / 0.2e1 * t92 + t120 / 0.2e1 + 0.3e1 / 0.2e1 * t104 + t96 / 0.2e1 - t108 / 0.2e1 + 0.2e1 * t125 - 0.2e1 * t127;
  t135 = sqrt(t93);
  t139 = log((-0.1e1 + t135) / t135);
  t141 = sqrt(t105);
  t145 = log((-0.1e1 + t141) / t141);
  
  double ret = 0.0;
  if (B_1 != 0.0)
    {
      ret += (t73 * CF + (0.67e2 / 0.9e1 - Pi2 + t78) * CA + (0.4e1 + 0.4e1 * t77) * beta0 - 0.10e2 / 0.9e1 * Nf) * B_1;
    }
  if (Bt_1 != 0.0)
    {
      ret += (t129 * CA + (-0.2e1 * t92 + 0.2e1 * t104 - 0.2e1 * t96 + 0.2e1 * t108 + 0.4e1 * t139 - 0.4e1 * t145 + 0.4e1 * (t141 - t135) / (0.1e1 + t135) / (0.1e1 + t141)) * beta0) * Bt_1;
    }
  ret *= AlphaS/TwoPi;
  return ret;
}
