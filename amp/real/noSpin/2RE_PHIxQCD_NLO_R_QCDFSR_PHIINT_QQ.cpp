
#include "AMP_HEADER.h"

double Eval_R_FSR_INT_QQ (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);



  double t1;
  double t10;
  double t101;
  double t105;
  double t11;
  double t115;
  double t116;
  double t117;
  double t12;
  double t120;
  double t121;
  double t126;
  double t127;
  double t13;
  double t134;
  double t14;
  double t140;
  double t143;
  double t145;
  double t154;
  double t155;
  double t158;
  double t161;
  double t164;
  double t167;
  double t173;
  double t177;
  double t189;
  double t193;
  double t197;
  double t2;
  double t20;
  double t208;
  double t209;
  double t219;
  double t220;
  double t222;
  double t230;
  double t233;
  double t234;
  double t236;
  double t239;
  double t241;
  double t242;
  double t244;
  double t247;
  double t248;
  double t250;
  double t26;
  double t27;
  double t272;
  double t288;
  double t289;
  double t3;
  double t305;
  double t31;
  double t311;
  double t32;
  double t35;
  double t37;
  double t4;
  double t40;
  double t41;
  double t45;
  double t46;
  double t5;
  double t50;
  double t53;
  double t54;
  double t56;
  double t60;
  double t63;
  double t70;
  double t71;
  double t75;
  double t76;
  double t8;
  double t80;
  double t84;
  double t85;
  double t86;
  double t88;
  double t90;
  double t93;
  double t97;
  t1 = sp(k1, p3);
  t2 = EPS_(k2, p1, p2, p3);
  t3 = t1 * t2;
  t4 = sp(p2, p3);
  t5 = t4 * At_fA_re;
  t8 = t4 * Bt_fH_re;
  t10 = sp(p1, p2);
  t11 = sp(p1, p3);
  t12 = t10 * t11;
  t13 = sp(k2, p3);
  t14 = t13 * t13;
  t20 = t1 * t1;
  t26 = t10 * t14;
  t27 = sp(k1, p1);
  t31 = sp(k1, p2);
  t32 = t31 * At_fH_re;
  t35 = t4 * At_fH_re;
  t37 = t4 * Bt_fA_re;
  t40 = t10 * t20;
  t41 = sp(k2, p1);
  t45 = sp(k2, p2);
  t46 = t45 * At_fH_re;
  t50 = -t12 * t14 * At_fH_re + 0.2e1 * t12 * t14 * Bt_fA_re - t12 * t20 * At_fH_re + 0.2e1 * t12 * t20 * Bt_fA_re - 0.2e1 * t26 * t27 * At_fH_re - 0.2e1 * t40 * t41 * At_fH_re - 0.2e1 * t26 * t32 - t26 * t35 + 0.2e1 * t26 * t37 - 0.2e1 * t3 * t5 + t3 * t8 - t35 * t40 - 0.2e1 * t40 * t46;
  t53 = t11 * t11;
  t54 = t53 * t13;
  t56 = t31 * Bt_fA_re;
  t60 = t45 * Bt_fA_re;
  t63 = t53 * t1;
  t70 = t11 * t13;
  t71 = t31 * t31;
  t75 = t11 * t1;
  t76 = t45 * t45;
  t80 = t27 * t27;
  t84 = t13 * t27;
  t85 = t4 * t4;
  t86 = t85 * At_fH_re;
  t88 = -0.2e1 * t13 * t35 * t80 - 0.2e1 * t70 * t71 * At_fH_re - 0.2e1 * t75 * t76 * At_fH_re + t32 * t54 - t32 * t63 + 0.2e1 * t37 * t40 - t46 * t54 + t46 * t63 + 0.2e1 * t54 * t56 + 0.2e1 * t54 * t60 + 0.2e1 * t56 * t63 + 0.2e1 * t60 * t63 + t84 * t86;
  t90 = t85 * Bt_fA_re;
  t93 = t13 * t41;
  t97 = t1 * t27;
  t101 = t41 * t41;
  t105 = t1 * t41;
  t115 = t10 * t13;
  t116 = EPS_(k1, k2, p2, p3);
  t117 = t116 * At_fA_re;
  t120 = EPS_(k1, k2, p1, p3);
  t121 = t120 * At_fA_re;
  t126 = -0.2e1 * t1 * t101 * t35 - 0.2e1 * t1 * t12 * At_fH_re - 0.2e1 * t12 * t13 * At_fH_re + t105 * t86 + 0.2e1 * t105 * t90 - 0.4e1 * t115 * t117 - 0.4e1 * t115 * t121 - 0.2e1 * t115 * t35 + 0.2e1 * t84 * t90 - t86 * t93 - t86 * t97 + 0.2e1 * t90 * t93 + 0.2e1 * t90 * t97;
  t127 = t10 * t1;
  t134 = EPS_(k1, p1, p2, p3);
  t140 = t2 * At_fA_re;
  t143 = t2 * Bt_fH_re;
  t145 = t11 * t134;
  t154 = t13 * t134;
  t155 = t27 * At_fA_re;
  t158 = t31 * At_fA_re;
  t161 = -0.2e1 * t1 * t145 * At_fA_re - t1 * t145 * Bt_fH_re + 0.2e1 * t134 * t70 * At_fA_re - t134 * t70 * Bt_fH_re + 0.4e1 * t117 * t127 + 0.4e1 * t121 * t127 - 0.2e1 * t127 * t35 - 0.2e1 * t140 * t70 + 0.2e1 * t140 * t75 - t143 * t70 - t143 * t75 + 0.4e1 * t154 * t155 - 0.4e1 * t154 * t158;
  t164 = t41 * At_fA_re;
  t167 = t45 * At_fA_re;
  t173 = t13 * t2;
  t177 = t134 * t1;
  t189 = t13 * t1;
  t193 = 0.4e1 * t12 * t189 * Bt_fA_re - 0.4e1 * t154 * t164 + 0.4e1 * t154 * t167 - 0.2e1 * t154 * t5 + t154 * t8 - 0.4e1 * t155 * t3 + 0.4e1 * t158 * t3 + 0.4e1 * t164 * t3 - 0.4e1 * t167 * t3 + 0.2e1 * t173 * t5 + t173 * t8 + 0.2e1 * t177 * t5 + t177 * t8;
  t197 = t1 * t31;
  t208 = t1 * t4;
  t209 = t208 * At_fH_re;
  t219 = t27 * t4;
  t220 = t219 * At_fH_re;
  t222 = t219 * Bt_fA_re;
  t230 = t31 * t45 * At_fH_re;
  t233 = t31 * t4;
  t234 = t233 * At_fH_re;
  t236 = t233 * Bt_fA_re;
  t239 = 0.2e1 * t1 * t115 * t45 * At_fH_re + 0.2e1 * t27 * t31 * t70 * At_fH_re - 0.2e1 * t31 * t41 * t70 * At_fH_re + 0.2e1 * t105 * t115 * At_fH_re + 0.2e1 * t115 * t197 * At_fH_re + 0.4e1 * t115 * t208 * Bt_fA_re + 0.2e1 * t115 * t97 * At_fH_re + 0.2e1 * t115 * t209 - t220 * t70 - 0.2e1 * t222 * t70 + 0.2e1 * t230 * t70 - t234 * t70 - 0.2e1 * t236 * t70;
  t241 = t41 * t4;
  t242 = t241 * At_fH_re;
  t244 = t241 * Bt_fA_re;
  t247 = t45 * t4;
  t248 = t247 * At_fH_re;
  t250 = t247 * Bt_fA_re;
  t272 = -0.2e1 * t27 * t45 * t75 * At_fH_re + 0.2e1 * t41 * t45 * t75 * At_fH_re + t220 * t75 - 0.2e1 * t222 * t75 + 0.2e1 * t230 * t75 + t234 * t75 - 0.2e1 * t236 * t75 + t242 * t70 - t242 * t75 - 0.2e1 * t244 * t70 - 0.2e1 * t244 * t75 + t248 * t70 - 0.2e1 * t250 * t70;
  t288 = sp(k1, k2);
  t289 = t288 * t10;
  t305 = 0.2e1 * t13 * t289 * t4 * At_fH_re + 0.2e1 * t12 * t189 * At_fH_re + 0.2e1 * t289 * t70 * At_fH_re + 0.2e1 * t289 * t75 * At_fH_re + 0.2e1 * t105 * t248 - 0.2e1 * t197 * t242 + 0.2e1 * t209 * t289 + 0.2e1 * t234 * t84 + 0.2e1 * t242 * t84 + 0.2e1 * t242 * t97 - t248 * t75 - 0.2e1 * t248 * t84 - 0.2e1 * t250 * t75;
  t311 = t10 * t10;
  return(-0.32e2 * PREF_R * (t50 + t88 + t126 + t161 + t193 + t239 + t272 + t305) / t1 / t311 / t13);
}