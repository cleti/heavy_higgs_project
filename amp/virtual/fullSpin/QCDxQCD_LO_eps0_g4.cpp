
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_GG (AMP_ARGS)
{


  AMP_DEFINITIONS;
  AP_REFS_B(ap);
  

  double t1;
  double t103;
  double t104;
  double t105;
  double t106;
  double t107;
  double t11;
  double t113;
  double t114;
  double t115;
  double t116;
  double t12;
  double t123;
  double t124;
  double t13;
  double t131;
  double t136;
  double t139;
  double t14;
  double t143;
  double t147;
  double t15;
  double t150;
  double t153;
  double t156;
  double t16;
  double t167;
  double t169;
  double t172;
  double t19;
  double t2;
  double t20;
  double t203;
  double t21;
  double t228;
  double t236;
  double t237;
  double t239;
  double t24;
  double t242;
  double t245;
  double t253;
  double t260;
  double t286;
  double t29;
  double t294;
  double t3;
  double t30;
  double t31;
  double t32;
  double t35;
  double t36;
  double t37;
  double t4;
  double t40;
  double t43;
  double t44;
  double t45;
  double t46;
  double t49;
  double t5;
  double t52;
  double t55;
  double t59;
  double t6;
  double t60;
  double t64;
  double t68;
  double t7;
  double t71;
  double t72;
  double t81;
  double t84;
  double t85;
  double t94;
  double t95;
  double t97;
  t1 = sp(s1, p1);
  t2 = sp(s2, p1);
  t3 = t1 * t2;
  t4 = beta * beta;
  t5 = t4 * t4;
  t6 = PREF_B_QCDxQCD_CFCA * t5;
  t7 = y * y;
  t11 = t1 * PREF_B_QCDxQCD_CA;
  t12 = sp(k1, s2);
  t13 = t4 * beta;
  t14 = t12 * t13;
  t15 = t7 * y;
  t16 = t14 * t15;
  t19 = t1 * PREF_B_QCDxQCD_CF;
  t20 = t12 * t5;
  t21 = t20 * t7;
  t24 = t1 * PREF_B_QCDxQCD_CFCA;
  t29 = t2 * PREF_B_QCDxQCD_CA;
  t30 = sp(s1, k2);
  t31 = t30 * t13;
  t32 = t31 * t15;
  t35 = t2 * PREF_B_QCDxQCD_CF;
  t36 = t30 * t5;
  t37 = t36 * t7;
  t40 = t2 * PREF_B_QCDxQCD_CFCA;
  t43 = sp(s1, s2);
  t44 = PREF_B_QCDxQCD_CA * t43;
  t45 = t5 * s;
  t46 = t45 * t7;
  t49 = PREF_B_QCDxQCD_CF * t43;
  t52 = PREF_B_QCDxQCD_CFCA * t43;
  t55 = PREF_B_QCDxQCD_CA * t4;
  t59 = -0.16e2 * t3 * t55 * t7 - 0.16e2 * t3 * t6 * t7 + 0.4e1 * t11 * t16 + 0.16e2 * t16 * t24 - 0.8e1 * t19 * t21 + 0.8e1 * t21 * t24 - 0.4e1 * t29 * t32 - 0.8e1 * t35 * t37 + 0.8e1 * t37 * t40 - 0.2e1 * t44 * t46 - 0.4e1 * t46 * t49 - 0.4e1 * t46 * t52;
  t60 = PREF_B_QCDxQCD_CF * t4;
  t64 = PREF_B_QCDxQCD_CFCA * t4;
  t68 = t14 * y;
  t71 = t12 * t4;
  t72 = t71 * t7;
  t81 = t31 * y;
  t84 = t30 * t4;
  t85 = t84 * t7;
  t94 = s * t4;
  t95 = t94 * t7;
  t97 = 0.16e2 * t3 * t60 * t7 - 0.48e2 * t3 * t64 * t7 - 0.4e1 * t11 * t68 + 0.8e1 * t11 * t72 - 0.8e1 * t19 * t72 - 0.16e2 * t24 * t68 + 0.24e2 * t24 * t72 + 0.4e1 * t29 * t81 + 0.8e1 * t29 * t85 - 0.8e1 * t35 * t85 + 0.16e2 * t40 * t81 + 0.24e2 * t40 * t85 - t44 * t95;
  t103 = t5 * t4;
  t104 = t103 * s;
  t105 = t7 * t7;
  t106 = t105 * t7;
  t107 = t104 * t106;
  t113 = t5 * beta;
  t114 = t12 * t113;
  t115 = t105 * y;
  t116 = t114 * t115;
  t123 = t30 * t113;
  t124 = t123 * t115;
  t131 = t104 * t105;
  t136 = -t107 * t44 + 0.4e1 * t107 * t49 - 0.4e1 * t107 * t52 - 0.4e1 * t11 * t116 + 0.16e2 * t116 * t19 - 0.16e2 * t116 * t24 + 0.4e1 * t124 * t29 - 0.16e2 * t124 * t35 + 0.16e2 * t124 * t40 + 0.2e1 * t131 * t44 - 0.8e1 * t131 * t49 + 0.4e1 * t49 * t95 - 0.4e1 * t52 * t95;
  t139 = PREF_B_QCDxQCD_CA * t5;
  t143 = PREF_B_QCDxQCD_CF * t5;
  t147 = t114 * t15;
  t150 = t20 * t105;
  t153 = t123 * t15;
  t156 = t36 * t105;
  t167 = t45 * t105;
  t169 = t104 * t7;
  t172 = 0.8e1 * t105 * t139 * t3 - 0.32e2 * t105 * t143 * t3 + 0.8e1 * t131 * t52 + 0.16e2 * t147 * t24 - 0.16e2 * t150 * t24 - 0.4e1 * t153 * t29 + 0.16e2 * t153 * t35 - 0.16e2 * t153 * t40 - 0.4e1 * t156 * t29 + 0.16e2 * t156 * t35 - 0.16e2 * t156 * t40 + t167 * t44 + 0.4e1 * t169 * t49;
  t203 = 0.32e2 * t105 * t3 * t6 + 0.16e2 * t143 * t3 * t7 - 0.4e1 * s * t143 - 0.4e1 * s * t6 + 0.4e1 * s * t60 + 0.4e1 * s * t64 + 0.4e1 * t11 * t147 - 0.4e1 * t11 * t150 - 0.16e2 * t147 * t19 + 0.16e2 * t150 * t19 + 0.4e1 * t167 * t52 - 0.4e1 * t169 * t52 - 0.16e2 * t32 * t40;
  t228 = s * t44 + 0.4e1 * s * t52 - PREF_B_QCDxQCD_CA * s + 0.4e1 * PREF_B_QCDxQCD_CF * s - 0.4e1 * t11 * t12 + 0.8e1 * t12 * t19 - 0.8e1 * t12 * t24 - 0.4e1 * t29 * t30 + 0.8e1 * t3 * PREF_B_QCDxQCD_CA - 0.16e2 * t3 * PREF_B_QCDxQCD_CF + 0.16e2 * t3 * PREF_B_QCDxQCD_CFCA + 0.8e1 * t30 * t35 - 0.8e1 * t30 * t40;
  t236 = PREF_B_QCDxQCD_CA * t103;
  t237 = s * t106;
  t239 = PREF_B_QCDxQCD_CF * t103;
  t242 = PREF_B_QCDxQCD_CFCA * t103;
  t245 = s * t105;
  t253 = s * t7;
  t260 = -t139 * t245 - 0.4e1 * t143 * t245 + t236 * t237 - 0.2e1 * t236 * t245 - 0.4e1 * t237 * t239 + 0.4e1 * t237 * t242 + 0.8e1 * t239 * t245 - 0.4e1 * t239 * t253 - 0.8e1 * t242 * t245 + 0.4e1 * t242 * t253 - 0.8e1 * t40 * t84 - 0.4e1 * t49 * t94 - 0.4e1 * t52 * t94;
  t286 = 0.2e1 * t139 * t253 + 0.4e1 * t143 * t253 - 0.8e1 * t19 * t71 - 0.8e1 * t24 * t71 + t253 * t55 + 0.4e1 * t253 * t6 - 0.4e1 * t253 * t60 - 0.4e1 * t253 * t64 + 0.16e2 * t3 * t60 + 0.16e2 * t3 * t64 - 0.8e1 * t35 * t84 + 0.4e1 * t45 * t49 + 0.4e1 * t45 * t52;
  t294 = pow(t4 * t7 - 0.1e1, 0.2e1);
  return(0.16e2 * (t59 + t97 + t136 + t172 + t203 + t228 + t260 + t286) / s / t294);
}
