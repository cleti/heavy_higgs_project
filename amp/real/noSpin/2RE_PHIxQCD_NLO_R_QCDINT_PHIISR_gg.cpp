
#include "AMP_HEADER.h"

double Eval_R_INT_ISR (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t100;
  double t103;
  double t107;
  double t108;
  double t109;
  double t11;
  double t113;
  double t114;
  double t115;
  double t116;
  double t117;
  double t118;
  double t119;
  double t12;
  double t120;
  double t127;
  double t13;
  double t136;
  double t137;
  double t139;
  double t14;
  double t140;
  double t143;
  double t149;
  double t150;
  double t153;
  double t157;
  double t158;
  double t159;
  double t160;
  double t166;
  double t167;
  double t17;
  double t178;
  double t179;
  double t18;
  double t180;
  double t181;
  double t182;
  double t183;
  double t186;
  double t189;
  double t190;
  double t191;
  double t192;
  double t193;
  double t194;
  double t195;
  double t196;
  double t197;
  double t198;
  double t199;
  double t2;
  double t200;
  double t201;
  double t202;
  double t207;
  double t208;
  double t209;
  double t210;
  double t212;
  double t22;
  double t227;
  double t23;
  double t232;
  double t238;
  double t241;
  double t249;
  double t25;
  double t250;
  double t251;
  double t253;
  double t256;
  double t257;
  double t258;
  double t259;
  double t26;
  double t27;
  double t274;
  double t275;
  double t28;
  double t284;
  double t285;
  double t29;
  double t298;
  double t30;
  double t302;
  double t303;
  double t31;
  double t316;
  double t32;
  double t321;
  double t322;
  double t327;
  double t328;
  double t336;
  double t341;
  double t347;
  double t35;
  double t354;
  double t355;
  double t360;
  double t363;
  double t367;
  double t368;
  double t375;
  double t38;
  double t382;
  double t383;
  double t384;
  double t385;
  double t390;
  double t396;
  double t397;
  double t4;
  double t40;
  double t41;
  double t411;
  double t42;
  double t422;
  double t44;
  double t440;
  double t459;
  double t47;
  double t495;
  double t500;
  double t51;
  double t511;
  double t513;
  double t516;
  double t52;
  double t520;
  double t523;
  double t525;
  double t526;
  double t528;
  double t529;
  double t53;
  double t530;
  double t531;
  double t533;
  double t538;
  double t539;
  double t545;
  double t546;
  double t547;
  double t549;
  double t55;
  double t550;
  double t554;
  double t57;
  double t571;
  double t58;
  double t585;
  double t588;
  double t6;
  double t61;
  double t617;
  double t62;
  double t630;
  double t634;
  double t637;
  double t638;
  double t641;
  double t662;
  double t69;
  double t7;
  double t717;
  double t719;
  double t738;
  double t750;
  double t755;
  double t758;
  double t772;
  double t78;
  double t8;
  double t83;
  double t843;
  double t845;
  double t846;
  double t847;
  double t849;
  double t86;
  double t9;
  double t93;
  double t94;
  double t97;
  t1 = sp(k1, p1);
  t2 = t1 * t1;
  t4 = sp(p2, p3);
  t6 = sp(k1, p2);
  t7 = 0.1e1 / t6;
  t8 = sp(k2, p1);
  t9 = 0.1e1 / t8;
  t11 = sp(p1, p2);
  t12 = t11 * t11;
  t13 = 0.1e1 / t12;
  t14 = t7 * t9 * t13;
  t17 = sp(k2, p2);
  t18 = 0.1e1 / t17;
  t22 = sp(p1, p3);
  t23 = 0.1e1 / t22;
  t25 = t4 * t4;
  t26 = t22 * t22;
  t27 = t25 + t26;
  t28 = t27 * t23;
  t29 = t4 + t22;
  t30 = 0.1e1 / t29;
  t31 = 0.1e1 / t11;
  t32 = t30 * t31;
  t35 = t25 * t13;
  t38 = 0.128e3 * t23 * t35 - 0.256e3 * t28 * t32 + 0.256e3 * t23;
  t40 = t4 * t22;
  t41 = 0.3e1 * t40;
  t42 = 0.2e1 * t26;
  t44 = (t25 - t41 - t42) * t23;
  t47 = sp(k1, p3);
  t51 = t4 - t22;
  t52 = t51 * t23;
  t53 = t30 * t11;
  t55 = 0.64e2 * t52 * t53;
  t57 = t29 * t29;
  t58 = 0.1e1 / t57;
  t61 = 0.3e1 * t22;
  t62 = t4 - t61;
  t69 = 0.128e3 * t40 * t13;
  t78 = 0.1e1 / t4;
  t83 = t26 * t13;
  t86 = 0.256e3 * t27 * t32 * t78 - 0.128e3 * t78 * t83 - 0.256e3 * t78;
  t93 = t30 * t78;
  t94 = t93 * t11;
  t97 = t23 * t58;
  t100 = t51 * t30;
  t103 = 0.256e3 * t27 * t51 * t78 * t97 + 0.256e3 * t100 * t31 - 0.256e3 * t52 * t94;
  t107 = 0.2e1 * t25;
  t108 = 0.3e1 * t26;
  t109 = t107 + t41 + t108;
  t113 = t25 * t4;
  t114 = 0.5e1 * t113;
  t115 = t22 * t25;
  t116 = 0.5e1 * t115;
  t117 = t26 * t4;
  t118 = 0.2e1 * t117;
  t119 = t26 * t22;
  t120 = 0.2e1 * t119;
  t127 = 0.256e3 * t40;
  t136 = 0.64e2 * (t26 + 0.6e1 * t40 + t25) * t23 * t94;
  t137 = t51 * t51;
  t139 = 0.128e3 * t137 * t58;
  t140 = 0.2e1 * t4;
  t143 = t23 * t31;
  t149 = t107 - t41 - t26;
  t150 = t149 * t23;
  t153 = t47 * t47;
  t157 = t93 * t12;
  t158 = t52 * t157;
  t159 = 0.64e2 * t158;
  t160 = 0.2e1 * t22;
  t166 = 0.5e1 * t25;
  t167 = 0.2e1 * t40;
  t178 = t30 * t12;
  t179 = 0.128e3 * t178;
  t180 = 0.4e1 * t25;
  t181 = 0.8e1 * t40;
  t182 = 0.4e1 * t26;
  t183 = 0.3e1 * t119;
  t186 = t58 * t11;
  t189 = 0.4e1 * t113;
  t190 = 0.4e1 * t115;
  t191 = 0.4e1 * t117;
  t192 = 0.4e1 * t119;
  t193 = t25 * t25;
  t194 = 0.3e1 * t193;
  t195 = t22 * t113;
  t196 = 0.6e1 * t195;
  t197 = t25 * t26;
  t198 = 0.18e2 * t197;
  t199 = t119 * t4;
  t200 = 0.10e2 * t199;
  t201 = t26 * t26;
  t202 = 0.3e1 * t201;
  t207 = 0.2e1 * t193;
  t208 = 0.7e1 * t195;
  t209 = 0.2e1 * t197;
  t210 = 0.2e1 * t201;
  t212 = (t189 + t190 - t191 - t192 + t207 + t208 - t209 - t199 - t210) * t23;
  t227 = t6 * t6;
  t232 = (t107 + t41 - t26) * t78;
  t238 = 0.3e1 * t4;
  t241 = t78 * t31;
  t249 = 0.64e2 * t51 * t78 * t53;
  t250 = 0.3e1 * t25;
  t251 = t250 + t41 + t42;
  t253 = t78 * t58;
  t256 = 0.2e1 * t113;
  t257 = 0.2e1 * t115;
  t258 = 0.5e1 * t117;
  t259 = 0.5e1 * t119;
  t274 = t25 + t41 - t42;
  t275 = t274 * t23;
  t284 = 0.5e1 * t193;
  t285 = 0.7e1 * t197;
  t298 = t78 * t11;
  t302 = 0.6e1 * t199;
  t303 = 0.5e1 * t201;
  t316 = 0.7e1 * t119;
  t321 = 0.10e2 * t195;
  t322 = 0.2e1 * t199;
  t327 = 0.5e1 * t195;
  t328 = 0.7e1 * t199;
  t336 = 0.7e1 * t113;
  t341 = 0.2e1 * t195;
  t347 = 0.5e1 * t199;
  t354 = t52 * t30;
  t355 = t153 * t47;
  t360 = 0.128e3 * t158;
  t363 = 0.256e3 * t100;
  t367 = 0.4e1 * t4;
  t368 = 0.4e1 * t22;
  t375 = t253 * t11;
  t382 = t12 * t11;
  t383 = t30 * t382;
  t384 = 0.128e3 * t383;
  t385 = 0.4e1 * t40;
  t390 = 0.6e1 * t115;
  t396 = 0.3e1 * t113;
  t397 = 0.3e1 * t117;
  t411 = t18 * t13;
  t422 = t238 - t22;
  t440 = 0.5e1 * t26;
  t459 = t189 + t190 - t191 - t192 + t207 + t195 + t209 - t328 - t210;
  t495 = 0.6e1 * t117;
  t500 = 0.3e1 * t115;
  t511 = 0.1e1 / t1;
  t513 = 0.128e3 * t2 * t1 * t4 * t14 + (-0.128e3 * t4 * t18 * t13 + t38 * t9 + (0.128e3 * t44 * t30 * t31 * t9 * t47 + (-0.64e2 * t23 * t25 * t30 * t31 * t62 + 0.128e3 * t4 * t51 * t58 + t55 + t69) * t9) * t7) * t2 + ((-t18 * t38 + t86 * t9) * t6 + (-0.128e3 * t18 * t32 * t44 + t103 * t9) * t47 + (t55 - 0.128e3 * t4 * t109 * t97 + 0.64e2 * (t114 + t116 + t118 - t120) * t23 * t32 + (-0.128e3 * t113 * t23 - t127) * t13) * t18 + (-t136 + t139 - 0.64e2 * (t140 + t61) * t29 * t143 + 0.256e3 * t35) * t9 + (-0.128e3 * t150 * t30 * t78 * t9 * t153 + (-t159 - 0.128e3 * (t4 + t160) * t51 * t97 * t11 + 0.64e2 * t4 * (t166 - t167 + t26) * t97 + 0.128e3 * t149 * t30 * t31) * t9 * t47 + (t179 + 0.32e2 * (-t180 - t181 - t182 + t113 - t116 - t117 - t183) * t23 * t186 - 0.32e2 * (-t189 - t190 - t191 - t192 + t194 + t196 + t198 + t200 + t202) * t23 * t58 - 0.32e2 * t212 * t32) * t9) * t7) * t1 + (-0.128e3 * t13 * t22 * t9 - t18 * t86) * t227 + ((0.128e3 * t232 * t32 * t9 - t103 * t18) * t47 + (-t136 + t139 - 0.64e2 * t29 * (t238 + t160) * t241 + 0.256e3 * t83) * t18 + (-t249 - 0.128e3 * t22 * t251 * t253 - 0.64e2 * (t256 - t257 - t258 - t259) * t78 * t32 + (-0.128e3 * t119 * t78 - t127) * t13) * t9) * t6 + (0.128e3 * t150 * t18 * t93 - 0.128e3 * t275 * t9 * t93) * t153 + ((-t159 + 0.128e3 * t251 * t23 * t186 - 0.64e2 * (t284 + t196 + t285 - 0.4e1 * t199 - t210) * t23 * t253 + 0.128e3 * (t113 - t115 + t191 + t120) * t23 * t32) * t18 + (t159 + 0.128e3 * t109 * t58 * t298 + 0.64e2 * (t207 + 0.4e1 * t195 - t285 - t302 - t303) * t23 * t253 + 0.128e3 * (t256 + t190 - t117 + t119) * t78 * t32) * t9) * t47 + (t179 + 0.32e2 * (t180 + t181 + t182 + t113 - t115 - t117 - t316) * t23 * t186 - 0.32e2 * (t189 + t190 + t191 + t192 + t194 + t321 + t209 - t322 - t303) * t23 * t58 + 0.32e2 * (t189 + t190 - t191 - t192 + t207 + t327 - t209 - t328 - 0.10e2 * t201) * t23 * t32) * t18 + (t179 - 0.32e2 * (-t180 - t181 - t182 + t336 + t115 + t117 - t119) * t58 * t298 + 0.32e2 * (-t189 - t190 - t191 - t192 + t284 + t341 - t209 - t200 - t202) * t78 * t58 - 0.32e2 * (t189 + t190 - t191 - t192 + 0.10e2 * t193 + t208 + t209 - t347 - t210) * t78 * t32) * t9 + (0.128e3 * t354 * t298 * t9 * t355 + (-0.256e3 * t186 * t28 + t360 - t363) * t9 * t153 + (-0.32e2 * (-t367 + t368 + t25 + t167 + t26) * t23 * t157 + 0.32e2 * (-t189 + t190 - t191 + t192 + t194 + t321 + t198 + 0.14e2 * t199 + t202) * t23 * t375 + 0.32e2 * t212 * t93) * t9 * t47 + (-t384 + 0.64e2 * (t367 + t368 + t166 + t385 + t108) * t58 * t12 - 0.64e2 * (t180 + t182 + t114 + t390 + 0.7e1 * t117 + t120) * t58 * t11 + 0.64e2 * (t180 + t181 + t182 + t396 + t257 + t397) * t30 - 0.64e2 * (t367 + t107 + t26) * t4 * t31) * t9) * t7 + (0.128e3 * t227 * t6 * t22 * t411 + (-0.128e3 * t232 * t30 * t31 * t18 * t47 + (0.64e2 * t26 * t31 * t422 * t93 - 0.128e3 * t22 * t51 * t58 - t249 + t69) * t18) * t227 + (0.128e3 * t275 * t30 * t78 * t18 * t153 + (t159 + 0.128e3 * (t140 + t22) * t51 * t375 + 0.64e2 * t22 * (t25 - t167 + t440) * t253 - 0.128e3 * t274 * t30 * t31) * t18 * t47 + (t179 - 0.32e2 * (t180 + t181 + t182 + t396 + t115 + t258 - t119) * t78 * t186 - 0.32e2 * (-t189 - t190 - t191 - t192 + t194 + t321 + t198 + t302 + t202) * t78 * t58 + 0.32e2 * t459 * t78 * t32) * t18) * t6 - 0.128e3 * t354 * t298 * t18 * t355 + (-0.256e3 * t27 * t298 * t58 - t360 + t363) * t18 * t153 + (-0.32e2 * (t367 - t368 + t25 + t167 + t26) * t23 * t157 + 0.32e2 * (t189 - t190 + t191 - t192 + t194 + 0.14e2 * t195 + t198 + t200 + t202) * t23 * t375 - 0.32e2 * t459 * t23 * t93) * t18 * t47 + (-t384 + 0.64e2 * (t367 + t368 + t250 + t385 + t440) * t58 * t12 - 0.64e2 * (t180 + t182 + t256 + 0.7e1 * t115 + t495 + t259) * t58 * t11 + 0.64e2 * (t180 + t181 + t182 + t500 + t118 + t183) * t30 - 0.64e2 * (t368 + t25 + t42) * t22 * t31) * t18) * t511;
  t516 = EPS_(k1, p1, p2, p3);
  t520 = t516 * t18;
  t523 = t516 * t78;
  t525 = t23 * t516 + t523;
  t526 = 0.512e3 * t525;
  t528 = t22 * t516;
  t529 = t528 * t78;
  t530 = t4 * t516;
  t531 = t530 * t23;
  t533 = -0.256e3 * t529 - 0.256e3 * t531;
  t538 = t29 * t23 * t78;
  t539 = t516 * t31;
  t545 = t516 * t51 * t23;
  t546 = t545 * t94;
  t547 = 0.128e3 * t546;
  t549 = (t107 + t40 + t26) * t78;
  t550 = t58 * t516;
  t554 = t30 * t516 * t31;
  t571 = t516 * t9;
  t585 = t253 * t516;
  t588 = 0.5e1 * t40;
  t617 = t23 * t78;
  t630 = 0.128e3 * t545 * t157;
  t634 = t523 * t11;
  t637 = 0.12e2 * t115;
  t638 = 0.12e2 * t117;
  t641 = 0.12e2 * t197;
  t662 = (t25 + t40 + t42) * t23;
  t717 = t51 * (t107 + t40 + t42);
  t719 = 0.64e2 * t717 * t143;
  t738 = t12 * t78;
  t750 = 0.64e2 * t717 * t241;
  t755 = 0.128e3 * t52 * t93 * t382;
  t758 = t253 * t12;
  t772 = 0.256e3 * t383;
  t843 = 0.96e2 * t525;
  t845 = 0.96e2 * t531;
  t846 = 0.96e2 * t529;
  t847 = 0.64e2 * t516;
  t849 = 0.32e2 * t530 + 0.32e2 * t528;
  return(At_fH_re * t513 * PREF_R_CA + At_fA_re * (-0.256e3 * t2 * t516 * t14 + (0.256e3 * t520 * t13 + (t13 * t533 + t31 * t526) * t9 + (0.256e3 * t538 * t539 * t9 * t47 + (-0.256e3 * t13 * t528 - 0.128e3 * t44 * t554 + 0.256e3 * t549 * t550 + t547) * t9) * t7) * t1 + ((-t13 * t533 - t31 * t526) * t18 - 0.256e3 * t571 * t13) * t6 + (-0.256e3 * t18 * t538 * t539 + 0.256e3 * t538 * t539 * t9) * t47 + (t547 + 0.256e3 * (t113 + t116 + t191 + t120) * t23 * t585 - 0.128e3 * (t166 + t588 + t42) * t23 * t554 + (0.256e3 * t23 * t25 * t516 + 0.512e3 * t528) * t13) * t18 + (t547 - 0.256e3 * (t256 + t190 + t258 + t119) * t23 * t585 + 0.128e3 * (t107 + t588 + t440) * t78 * t554 + (-0.256e3 * t26 * t516 * t78 - 0.512e3 * t530) * t13) * t9 + (-0.256e3 * t617 * t571 * t153 + (-0.256e3 * t546 + 0.256e3 * (t256 - t115 - t118 - t183) * t23 * t585 + 0.256e3 * t549 * t554) * t9 * t47 + (-t630 + 0.64e2 * (t336 + t500 + t117 - t183) * t23 * t58 * t634 - 0.64e2 * (t189 + t637 + t638 + t192 + 0.9e1 * t193 + 0.12e2 * t195 + t641 - t201) * t23 * t585 - 0.64e2 * (t166 - t26) * t30 * t539) * t9) * t7 + (0.256e3 * t227 * t516 * t411 + (-0.256e3 * t538 * t539 * t18 * t47 + (0.256e3 * t13 * t530 - 0.128e3 * t232 * t554 - 0.256e3 * t550 * t662 + t547) * t18) * t6 + 0.256e3 * t617 * t520 * t153 + (-0.256e3 * t546 + 0.256e3 * (t396 + t257 + t117 - t120) * t23 * t585 - 0.256e3 * t662 * t554) * t18 * t47 + (-t630 + 0.64e2 * (t396 - t115 - t397 - t316) * t23 * t58 * t634 - 0.64e2 * (-t189 - t637 - t638 - t192 + t193 - t641 - 0.12e2 * t199 - 0.9e1 * t201) * t23 * t585 - 0.64e2 * (t25 - t440) * t30 * t539) * t18) * t511) * PREF_R_CA + Bt_fA_re * ((-0.128e3 * t12 * t23 + 0.64e2 * (t250 + t26) * t23 * t53 - 0.64e2 * (t396 - t115 + t117 + t119) * t23 * t30 + t719) * t9 * t7 * t1 + (0.128e3 * t62 * t23 * t178 - 0.64e2 * (t396 - t116 + t117 - t316) * t23 * t186 + 0.64e2 * (t396 - t116 + t117 - t183) * t23 * t30 - t719) * t18 + (-0.128e3 * t422 * t30 * t738 + 0.64e2 * (t336 - t115 + t258 - t183) * t58 * t298 - 0.64e2 * (t396 - t115 + t258 - t183) * t78 * t30 + t750) * t9 + ((t755 - 0.64e2 * (t396 - t116 - t397 - t183) * t23 * t758 + 0.64e2 * (t194 - t341 - t302 - t202) * t23 * t375 - 0.64e2 * (t207 + t195 - t347 - t210) * t23 * t93) * t9 * t47 + (t772 - 0.128e3 * (t166 + t385 + t108) * t58 * t12 + 0.128e3 * (t114 + t390 + 0.9e1 * t117 + t192) * t58 * t11 - 0.128e3 * (t396 + t257 + t258 + t120) * t30 + 0.128e3 * (t107 + t26) * t4 * t31) * t9) * t7 + ((-0.128e3 * t738 + 0.64e2 * (t25 + t108) * t78 * t53 - 0.64e2 * (t113 + t115 - t117 + t183) * t78 * t30 - t750) * t18 * t6 + (-t755 + 0.64e2 * (t396 + t500 + t258 - t183) * t23 * t758 - 0.64e2 * (t194 + t196 + t322 - t202) * t23 * t375 + 0.64e2 * (t207 + t327 - t199 - t210) * t23 * t93) * t18 * t47 + (t772 - 0.128e3 * (t250 + t385 + t440) * t58 * t12 + 0.128e3 * (t189 + 0.9e1 * t115 + t495 + t259) * t58 * t11 - 0.128e3 * (t256 + t116 + t118 + t183) * t30 + 0.128e3 * (t25 + t42) * t22 * t31) * t18) * t511) * PREF_R_CA + Bt_fH_re * ((t11 * t843 + t31 * t849 - t845 - t846 - t847) * t9 * t7 + (-t11 * t843 - t31 * t849 + t845 + t846 + t847) * t18 * t511) * PREF_R_CA);
}
