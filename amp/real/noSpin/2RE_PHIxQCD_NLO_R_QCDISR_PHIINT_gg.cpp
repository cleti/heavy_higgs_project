
#include "AMP_HEADER.h"

double Eval_R_ISR_INT (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);



  double t1;
  double t103;
  double t11;
  double t115;
  double t116;
  double t118;
  double t12;
  double t125;
  double t126;
  double t13;
  double t137;
  double t138;
  double t14;
  double t142;
  double t143;
  double t144;
  double t146;
  double t147;
  double t152;
  double t155;
  double t157;
  double t158;
  double t16;
  double t166;
  double t167;
  double t17;
  double t173;
  double t177;
  double t2;
  double t200;
  double t203;
  double t206;
  double t21;
  double t22;
  double t225;
  double t226;
  double t227;
  double t229;
  double t23;
  double t232;
  double t244;
  double t245;
  double t246;
  double t25;
  double t250;
  double t26;
  double t264;
  double t27;
  double t28;
  double t29;
  double t3;
  double t314;
  double t316;
  double t319;
  double t32;
  double t320;
  double t321;
  double t325;
  double t328;
  double t33;
  double t331;
  double t338;
  double t339;
  double t343;
  double t348;
  double t35;
  double t351;
  double t355;
  double t359;
  double t36;
  double t365;
  double t368;
  double t39;
  double t397;
  double t4;
  double t40;
  double t404;
  double t407;
  double t414;
  double t45;
  double t5;
  double t57;
  double t6;
  double t63;
  double t69;
  double t7;
  double t70;
  double t71;
  double t72;
  double t74;
  double t75;
  double t76;
  double t8;
  double t81;
  double t83;
  double t84;
  double t86;
  double t9;
  double t94;
  double t96;
  t1 = sp(p2, p3);
  t2 = sp(p1, p3);
  t3 = t1 - t2;
  t4 = t1 * t3;
  t5 = t1 + t2;
  t6 = t5 * t5;
  t7 = 0.1e1 / t6;
  t8 = 0.1e1 / t2;
  t9 = t7 * t8;
  t11 = 0.2e1 * t1;
  t12 = t11 + t2;
  t13 = t1 * t12;
  t14 = 0.1e1 / t5;
  t16 = sp(p1, p2);
  t17 = 0.1e1 / t16;
  t21 = -0.256e3 * t13 * t14 * t17 * t8 - 0.256e3 * t4 * t9;
  t22 = sp(k2, p2);
  t23 = 0.1e1 / t22;
  t25 = t1 * t14;
  t26 = t25 * t17;
  t27 = t16 * t16;
  t28 = 0.1e1 / t27;
  t29 = t1 * t28;
  t32 = sp(k2, p1);
  t33 = 0.1e1 / t32;
  t35 = sp(k1, p2);
  t36 = 0.1e1 / t35;
  t39 = sp(k1, p1);
  t40 = t39 * t39;
  t45 = 0.256e3 * (t2 + 0.7e1 * t1) * t7;
  t57 = 0.256e3 * (0.7e1 * t2 + t1) * t7;
  t63 = t2 * t28;
  t69 = t3 * t7;
  t70 = t8 * t16;
  t71 = t69 * t70;
  t72 = t1 * t1;
  t74 = t1 * t2;
  t75 = 0.2e1 * t74;
  t76 = t2 * t2;
  t81 = 0.256e3 * t71 + 0.256e3 * (0.3e1 * t72 + t75 + t76) * t7 * t8 + 0.256e3 * t26;
  t83 = 0.3e1 * t1;
  t84 = t2 + t83;
  t86 = 0.2e1 * t2;
  t94 = sp(k1, p3);
  t96 = t12 * t14;
  t103 = t14 * t17;
  t115 = t74 * t103;
  t116 = 0.256e3 * t115;
  t118 = 0.512e3 * t74 * t28;
  t125 = 0.4e1 * t74;
  t126 = t72 + t125 + t76;
  t137 = t2 * t14;
  t138 = t137 * t17;
  t142 = t2 * t3;
  t143 = 0.1e1 / t1;
  t144 = t7 * t143;
  t146 = t1 + t86;
  t147 = t2 * t146;
  t152 = -0.256e3 * t14 * t143 * t147 * t17 + 0.256e3 * t142 * t144;
  t155 = t35 * t35;
  t157 = 0.3e1 * t2;
  t158 = t1 + t157;
  t166 = t143 * t16;
  t167 = t69 * t166;
  t173 = -0.256e3 * t167 + 0.256e3 * (t72 + t75 + 0.3e1 * t76) * t7 * t143 + 0.256e3 * t138;
  t177 = t146 * t14;
  t200 = -0.256e3 * t71 - 0.256e3 * t14;
  t203 = 0.256e3 * t167 - 0.256e3 * t14;
  t206 = t94 * t94;
  t225 = t27 * t1 * t7;
  t226 = 0.512e3 * t225;
  t227 = 0.2e1 + t1 + t157;
  t229 = t7 * t16;
  t232 = 0.2e1 * t76;
  t244 = t27 * t2 * t7;
  t245 = 0.512e3 * t244;
  t246 = 0.2e1 + t83 + t2;
  t250 = 0.2e1 * t72;
  t264 = t126 * t7;
  t314 = 0.1e1 / t39;
  t316 = (t21 * t23 + (0.512e3 * t26 - 0.512e3 * t29) * t33 + t21 * t36) * t40 + (((t45 - 0.256e3 * (0.4e1 * t1 + t2) * t14 * t17 + 0.512e3 * t29) * t23 + (t57 - 0.256e3 * (t1 + 0.4e1 * t2) * t14 * t17 + 0.512e3 * t63) * t33) * t35 + (t81 * t23 + (-0.256e3 * t84 * t7 + 0.256e3 * (t83 + t86) * t14 * t17) * t33) * t94 + t45 - 0.256e3 * t96 * t17 + (-0.1024e4 * t1 * t16 * t7 - 0.1024e4 * t103 * t13 + 0.512e3 * t28 * t72 + 0.1024e4 * t25) * t23 + (0.256e3 * t2 * (0.5e1 * t1 - t2) * t7 + t116 - t118) * t33 + (0.256e3 * t1 * t126 * t9 - 0.256e3 * t14 * t17 * t72 + 0.256e3 * t16 * t4 * t9 + t81 * t94) * t36) * t39 + ((0.512e3 * t138 - 0.512e3 * t63) * t23 + t152 * t33) * t155 + (((-0.256e3 * t158 * t7 + 0.256e3 * (t11 + t157) * t14 * t17) * t23 + t173 * t33) * t94 + t57 - 0.256e3 * t177 * t17 + (-0.256e3 * t1 * (t1 - 0.5e1 * t2) * t7 + t116 - t118) * t23 + (-0.1024e4 * t16 * t2 * t7 - 0.1024e4 * t103 * t147 + 0.512e3 * t28 * t76 + 0.1024e4 * t137) * t33) * t35 + (t200 * t23 + t203 * t33) * t206 + (-0.1024e4 * t14 + 0.1280e4 * t17 + (0.512e3 * t4 * t7 + 0.512e3 * t115) * t23 + (-0.512e3 * t142 * t7 + 0.512e3 * t115) * t33) * t94 - 0.1024e4 * t14 * t16 + 0.1280e4 - 0.512e3 * t5 * t17 + (t226 - 0.512e3 * t1 * t227 * t229 + 0.512e3 * t1 * (t1 + t2 + t72 + t125 + t232) * t7 - 0.256e3 * t1 * (t11 + t86 + t72 + t125 + t232) * t103) * t23 + (t245 - 0.512e3 * t2 * t246 * t229 + 0.512e3 * t2 * (t1 + t2 + t250 + t125 + t76) * t7 - 0.256e3 * t2 * (t11 + t86 + t250 + t125 + t76) * t103) * t33 + (t200 * t206 + (-0.256e3 * t27 * t69 * t8 - 0.256e3 * t264 * t70 + 0.256e3 * t25) * t94 + t226 - 0.256e3 * t1 * (0.4e1 + t83 + t2) * t229 + 0.256e3 * t1 * t246 * t14 - 0.256e3 * (0.2e1 + t1) * t1 * t17) * t36 + (t152 * t155 + (0.256e3 * t126 * t144 * t2 - 0.256e3 * t14 * t17 * t76 - 0.256e3 * t142 * t144 * t16 + t173 * t94) * t35 + t203 * t206 + (0.256e3 * t143 * t27 * t69 - 0.256e3 * t166 * t264 + 0.256e3 * t137) * t94 + t245 - 0.256e3 * t2 * (0.4e1 + t1 + t157) * t229 + 0.256e3 * t2 * t227 * t14 - 0.256e3 * (0.2e1 + t2) * t2 * t17) * t314;
  t319 = EPS_(k1, p1, p2, p3);
  t320 = t8 * t319;
  t321 = t69 * t320;
  t325 = -0.512e3 * t17 * t320 * t96 - 0.512e3 * t321;
  t328 = t14 * t319 * t17;
  t331 = 0.1024e4 * t28 * t319 - 0.1024e4 * t328;
  t338 = t143 * t319;
  t339 = t69 * t338;
  t343 = 0.512e3 * t17 * t177 * t338 - 0.512e3 * t339;
  t348 = 0.512e3 * t321 + 0.512e3 * t328;
  t351 = 0.512e3 * t339 - 0.512e3 * t328;
  t355 = t1 * t319;
  t359 = t25 * t319 * t17;
  t365 = t2 * t319;
  t368 = t365 * t103;
  t397 = t1 * t84;
  t404 = -0.512e3 * t14 * t397 + 0.512e3 * t17 * t72 + 0.512e3 * t229 * t397 - 0.1024e4 * t225;
  t407 = t2 * t158;
  t414 = -0.512e3 * t14 * t407 + 0.512e3 * t17 * t76 + 0.512e3 * t229 * t407 - 0.1024e4 * t244;
  return(At_fH_re * t316 * PREF_R_CA + At_fA_re * ((t23 * t325 + t325 * t36 + t33 * t331) * t39 + (-t23 * t331 + t33 * t343) * t35 + (t23 * t348 + t33 * t351) * t94 + (-0.1024e4 * t28 * t355 - 0.2048e4 * t355 * t7 + 0.1024e4 * t359) * t23 + (0.1024e4 * t28 * t365 + 0.2048e4 * t365 * t7 - 0.1024e4 * t368) * t33 + (0.512e3 * t16 * t320 * t69 + 0.512e3 * t264 * t320 + t348 * t94 - 0.512e3 * t359) * t36 + (0.512e3 * t16 * t338 * t69 - 0.512e3 * t264 * t338 + t343 * t35 + t351 * t94 + 0.512e3 * t368) * t314) * PREF_R_CA + Bt_fA_re * (t23 * t404 + t314 * t414 + t33 * t414 + t36 * t404) * PREF_R_CA);
}