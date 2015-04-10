
#include "AMP_HEADER.h"

double Eval_V_B1 (AMP_ARGS)
{

AMP_DEFINITIONS


  c_double& cg  = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg1 = I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1;
  c_double& cg3 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0;
  c_double& cg5 = I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1;

  
  double t1;
  double t10;
  double t100;
  double t101;
  double t108;
  double t11;
  double t113;
  double t12;
  double t120;
  double t123;
  double t124;
  double t125;
  double t128;
  double t13;
  double t136;
  double t137;
  double t14;
  double t147;
  double t149;
  double t150;
  double t151;
  double t153;
  double t154;
  double t161;
  double t162;
  double t163;
  double t165;
  double t167;
  double t169;
  double t170;
  double t178;
  double t18;
  double t180;
  double t181;
  double t184;
  double t185;
  double t186;
  double t190;
  double t191;
  double t192;
  double t199;
  double t2;
  double t20;
  double t202;
  double t203;
  double t208;
  double t209;
  double t21;
  double t211;
  double t213;
  double t218;
  double t221;
  double t222;
  double t223;
  double t224;
  double t226;
  double t227;
  double t228;
  double t230;
  double t233;
  double t234;
  double t242;
  double t244;
  double t246;
  double t247;
  double t248;
  double t25;
  double t250;
  double t251;
  double t252;
  double t253;
  double t255;
  double t258;
  double t259;
  double t260;
  double t261;
  double t269;
  double t27;
  double t270;
  double t271;
  double t276;
  double t277;
  double t28;
  double t280;
  double t283;
  double t291;
  double t292;
  double t294;
  double t295;
  double t297;
  double t3;
  double t30;
  double t300;
  double t301;
  double t302;
  double t310;
  double t315;
  double t316;
  double t318;
  double t32;
  double t321;
  double t325;
  double t326;
  double t33;
  double t333;
  double t334;
  double t336;
  double t338;
  double t341;
  double t342;
  double t349;
  double t350;
  double t355;
  double t356;
  double t363;
  double t364;
  double t366;
  double t369;
  double t376;
  double t380;
  double t393;
  double t397;
  double t4;
  double t40;
  double t402;
  double t404;
  double t405;
  double t407;
  double t410;
  double t413;
  double t417;
  double t419;
  double t42;
  double t423;
  double t426;
  double t43;
  double t44;
  double t46;
  double t48;
  double t482;
  double t49;
  double t490;
  double t494;
  double t497;
  double t5;
  double t500;
  double t509;
  double t51;
  double t517;
  double t526;
  double t528;
  double t53;
  double t54;
  double t544;
  double t55;
  double t556;
  double t559;
  double t56;
  double t561;
  double t566;
  double t57;
  double t571;
  double t578;
  double t58;
  double t59;
  double t590;
  double t591;
  double t597;
  double t6;
  double t60;
  double t607;
  double t61;
  double t614;
  double t631;
  double t633;
  double t636;
  double t641;
  double t644;
  double t651;
  double t655;
  double t659;
  double t663;
  double t664;
  double t669;
  double t678;
  double t68;
  double t69;
  double t739;
  double t742;
  double t746;
  double t75;
  double t751;
  double t76;
  double t77;
  double t773;
  double t776;
  double t78;
  double t787;
  double t79;
  double t795;
  double t8;
  double t80;
  double t801;
  double t802;
  double t807;
  double t81;
  double t82;
  double t84;
  double t848;
  double t854;
  double t866;
  double t891;
  double t91;
  double t94;
  double t96;
  double t97;
  t1 = beta_y + 0.1e1;
  t2 = 0.1e1 / t1;
  t3 = beta_y - 0.1e1;
  t4 = 0.1e1 / t3;
  t5 = t2 * t4;
  t6 = sp(s1, k2);
  t8 = 0.1e1 / s;
  t10 = beta_y * t8 * t4;
  t11 = -s + 0.2e1;
  t12 = t2 * t11;
  t13 = sp(s1, p1);
  t14 = t12 * t13;
  t18 = sp(k1, s2);
  t20 = sp(s2, p1);
  t21 = t6 * t20;
  t25 = beta * beta;
  t27 = y * y;
  t28 = t25 * s * t27;
  t30 = 0.2e1 * t25 * t27;
  t32 = (-t28 + t30 + s - 0.4e1) * t2;
  t33 = sp(s1, s2);
  t40 = RE(I2_MT2_0_MT2_MU2_0);
  t42 = beta_y * s;
  t43 = 0.6e1 * t42;
  t44 = 0.5e1 * s;
  t46 = (t28 - t43 + t44 - 0.8e1) * t4;
  t48 = 0.1e1 / (t42 - s + 0.2e1);
  t49 = t48 * t6;
  t51 = s * s;
  t53 = t25 * t51 * t27;
  t54 = beta_y * t51;
  t55 = 0.2e1 * t54;
  t56 = 0.4e1 * t42;
  t57 = 0.3e1 * t51;
  t58 = 0.12e2 * s;
  t59 = t53 + t55 - t56 - t57 + t58 - 0.16e2;
  t60 = t59 * t8;
  t61 = t4 * t48;
  t68 = t48 * t20;
  t69 = t68 * t6;
  t75 = t25 * beta * t51 * t27 * y;
  t76 = 0.3e1 * t53;
  t77 = 0.8e1 * t28;
  t78 = 0.9e1 * t54;
  t79 = 0.36e2 * t42;
  t80 = 0.5e1 * t51;
  t81 = 0.16e2 * beta_y;
  t82 = 0.28e2 * s;
  t84 = (t75 + t76 - t77 - t78 + t79 + t80 - t81 - t82 + 0.32e2) * t4;
  t91 = RE(I2_T11_0_MT2_MU2_0);
  t94 = (t28 + t43 + t44 - 0.8e1) * t2;
  t96 = 0.1e1 / (t42 + s - 0.2e1);
  t97 = t96 * t6;
  t100 = (t53 - t55 + t56 - t57 + t58 - 0.16e2) * t2;
  t101 = t8 * t96;
  t108 = t96 * t20;
  t113 = (t75 - t76 + t77 - t78 + t79 - t80 - t81 + t82 - 0.32e2) * t2;
  t120 = RE(I2_T12_0_MT2_MU2_0);
  t123 = (t28 - s + 0.2e1) * s;
  t124 = t96 * t48;
  t125 = t124 * t6;
  t128 = t124 * t13;
  t136 = t28 - s + 0.4e1;
  t137 = t11 * t136;
  t147 = RE(I1_MT2_MU2_0);
  t149 = t3 * s;
  t150 = 0.1e1 / t136;
  t151 = t150 * t6;
  t153 = t42 - s + 0.4e1;
  t154 = t153 * t150;
  t161 = s * t33;
  t162 = 0.32e2 * t161;
  t163 = 0.32e2 * s;
  t165 = RE(I3_MT2_0_T11_0_MT2_MT2_MU2_0);
  t167 = t1 * s;
  t169 = t42 + s - 0.4e1;
  t170 = t169 * t150;
  t178 = RE(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
  t180 = -s + 0.4e1;
  t181 = s * t180;
  t184 = t136 * t136;
  t185 = 0.1e1 / t184;
  t186 = (0.3e1 * t28 - s + 0.4e1) * t185;
  t190 = t180 * t180;
  t191 = t190 * t185;
  t192 = t191 * t13;
  t199 = t185 * t6 * t20;
  t202 = t190 * s;
  t203 = t150 * t33;
  t208 = (0.64e2 * t181 * t186 * t6 - 0.128e3 * t192 * t42) * t18 + 0.128e3 * t42 * t190 * t199 + 0.32e2 * t202 * t203 - 0.32e2 * t202 * t150;
  t209 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t211 = s * t150;
  t213 = t150 * t13;
  t218 = t151 * t20;
  t221 = 0.64e2 * t161;
  t222 = 0.64e2 * s;
  t223 = (0.128e3 * t211 * t6 + 0.128e3 * t213 * t42) * t18 - 0.128e3 * t42 * t218 + t221 - t222;
  t224 = RE(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t226 = 0.4e1 * t28;
  t227 = 0.2e1 * t42;
  t228 = 0.6e1 * s;
  t230 = (-t53 + t226 - t227 + t51 - t228 + 0.8e1) * s;
  t233 = t180 * t3;
  t234 = t211 * t13;
  t242 = t11 * t180;
  t244 = 0.8e1 * t242 * t161;
  t246 = 0.8e1 * t242 * s;
  t247 = (0.16e2 * t151 * t230 - 0.32e2 * t233 * t234) * t18 + 0.32e2 * t233 * s * t218 + t244 - t246;
  t248 = RE(cg);
  t250 = t51 * t180;
  t251 = t3 * t3;
  t252 = t250 * t251;
  t253 = 0.2e1 * t28;
  t255 = (t253 - t42 - s + 0.4e1) * t185;
  t258 = t51 * t190;
  t259 = t251 * t3;
  t260 = t259 * t185;
  t261 = t260 * t13;
  t269 = t190 * t251;
  t270 = t51 * t150;
  t271 = t270 * t33;
  t276 = (0.16e2 * t252 * t255 * t6 - 0.16e2 * t258 * t261) * t18 + 0.16e2 * t258 * t259 * t199 + 0.8e1 * t269 * t271 - 0.8e1 * t269 * t270;
  t277 = RE(cg1);
  t280 = (-t53 + t226 + t227 + t51 - t228 + 0.8e1) * s;
  t283 = t180 * t1;
  t291 = (0.16e2 * t151 * t280 - 0.32e2 * t234 * t283) * t18 + 0.32e2 * t283 * s * t218 + t244 - t246;
  t292 = RE(cg3);
  t294 = t1 * t1;
  t295 = t250 * t294;
  t297 = (t253 + t42 - s + 0.4e1) * t185;
  t300 = t294 * t1;
  t301 = t300 * t185;
  t302 = t301 * t13;
  t310 = t190 * t294;
  t315 = (0.16e2 * t295 * t297 * t6 - 0.16e2 * t258 * t302) * t18 + 0.16e2 * t258 * t300 * t199 + 0.8e1 * t310 * t271 - 0.8e1 * t310 * t270;
  t316 = RE(cg5);
  t318 = EPS_(k1, k2, p1, s1);
  t321 = EPS_(k1, k2, p1, s2);
  t325 = 0.128e3 * t191 * t318 * t42 + 0.128e3 * t191 * t321 * t42;
  t326 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t333 = -0.128e3 * t150 * t318 * t42 - 0.128e3 * t150 * t321 * t42;
  t334 = IM(I3_S12_0_0_MT2_MT2_MT2_MU2_0);
  t336 = t211 * t318;
  t338 = t211 * t321;
  t341 = 0.32e2 * t233 * t336 + 0.32e2 * t233 * t338;
  t342 = IM(cg);
  t349 = 0.16e2 * t258 * t260 * t318 + 0.16e2 * t258 * t260 * t321;
  t350 = IM(cg1);
  t355 = 0.32e2 * t283 * t336 + 0.32e2 * t283 * t338;
  t356 = IM(cg3);
  t363 = 0.16e2 * t258 * t301 * t318 + 0.16e2 * t258 * t301 * t321;
  t364 = IM(cg5);
  t366 = ((-0.128e3 * t10 * t14 + 0.128e3 * t5 * t6) * t18 + 0.128e3 * t10 * t12 * t21 - 0.64e2 * t32 * t4 * t33 + 0.64e2 * t32 * t4) * t40 + ((-0.16e2 * t13 * t60 * t61 + 0.16e2 * t46 * t49) * t18 + 0.16e2 * t60 * t4 * t69 - 0.8e1 * t84 * t48 * t33 + 0.8e1 * t84 * t48) * t91 + ((0.16e2 * t100 * t101 * t13 + 0.16e2 * t94 * t97) * t18 - 0.16e2 * t100 * t8 * t108 * t6 + 0.8e1 * t113 * t96 * t33 - 0.8e1 * t113 * t96) * t120 + ((-0.32e2 * t123 * t125 + 0.64e2 * t128 * t42) * t18 - 0.64e2 * t42 * t96 * t69 - 0.16e2 * t137 * s * t124 * t33 + 0.16e2 * t137 * s * t96 * t48) * t147 + ((-0.64e2 * t13 * t154 + 0.64e2 * t149 * t151) * t18 + 0.64e2 * t154 * t21 - t162 + t163) * t165 + ((-0.64e2 * t13 * t170 - 0.64e2 * t151 * t167) * t18 + 0.64e2 * t170 * t21 - t162 + t163) * t178 + t208 * t209 + t223 * t224 + t247 * t248 + t276 * t277 + t291 * t292 + t315 * t316 + t325 * t326 + t333 * t334 + t341 * t342 + t349 * t350 + t355 * t356 + t363 * t364;
  t369 = EPS_(k1, k2, s1, s2);
  t376 = t150 * t369;
  t380 = t11 * s;
  t393 = t150 * t18;
  t397 = 0.128e3 * t151 * t181 + 0.128e3 * t181 * t393;
  t402 = 0.32e2 * t18 * t380 + 0.32e2 * t380 * t6;
  t404 = t251 * t180;
  t405 = t270 * t18;
  t407 = t270 * t6;
  t410 = 0.32e2 * t404 * t405 + 0.32e2 * t404 * t407;
  t413 = t180 * t294;
  t417 = 0.32e2 * t405 * t413 + 0.32e2 * t407 * t413;
  t419 = -0.128e3 * t181 * t209 * t376 - 0.32e2 * t248 * t369 * t380 - 0.32e2 * t252 * t277 * t376 - 0.32e2 * t292 * t369 * t380 - 0.32e2 * t295 * t316 * t376 - 0.64e2 * t120 * t369 + t326 * t397 + t342 * t402 + t350 * t410 + t356 * t402 + t364 * t417 + 0.128e3 * t369 * t40 - 0.64e2 * t369 * t91;
  t423 = EPS_(k1, p1, s1, s2);
  t426 = EPS_(k2, p1, s1, s2);
  t482 = 0.64e2 * t181 * t186 * t369 - 0.128e3 * t191 * t42 * t423 - 0.128e3 * t191 * t42 * t426;
  t490 = 0.128e3 * t150 * t42 * t423 + 0.128e3 * t150 * t42 * t426 + 0.128e3 * t211 * t369;
  t494 = t211 * t423;
  t497 = t211 * t426;
  t500 = 0.16e2 * t230 * t376 - 0.32e2 * t233 * t494 - 0.32e2 * t233 * t497;
  t509 = 0.16e2 * t252 * t255 * t369 - 0.16e2 * t258 * t260 * t423 - 0.16e2 * t258 * t260 * t426;
  t517 = 0.16e2 * t280 * t376 - 0.32e2 * t283 * t494 - 0.32e2 * t283 * t497;
  t526 = -0.16e2 * t258 * t301 * t423 - 0.16e2 * t258 * t301 * t426 + 0.16e2 * t295 * t297 * t369;
  t528 = 0.2e1 * s;
  t544 = 0.32e2 * t181 * (t53 + t253 + t54 - t56 - t528 + 0.8e1) * t185 * t18 + 0.32e2 * t181 * (t53 + t253 - t54 + t56 - t528 + 0.8e1) * t185 * t6 - 0.64e2 * t54 * t192 + 0.64e2 * t54 * t191 * t20;
  t556 = t150 * t20;
  t559 = 0.32e2 * (t28 + t42 + 0.4e1) * s * t393 + 0.32e2 * (t28 - t42 + 0.4e1) * s * t151 + 0.64e2 * t54 * t213 - 0.64e2 * t54 * t556;
  t561 = 0.8e1 * s;
  t566 = 0.16e2 * s;
  t571 = t3 * t150;
  t578 = 0.8e1 * (-t53 + t226 - t56 + t51 - t561 + 0.16e2) * s * t393 + 0.8e1 * (-t53 + t226 - t55 + t56 + t57 - t566 + 0.16e2) * s * t151 - 0.16e2 * t250 * t571 * t13 + 0.16e2 * t250 * t571 * t20;
  t590 = t51 * s;
  t591 = t590 * t190;
  t597 = 0.4e1 * t252 * (t53 + t226 - t56 - t51 + 0.16e2) * t185 * t18 + 0.4e1 * t252 * (t53 + t226 - t55 + t56 + t51 - t561 + 0.16e2) * t185 * t6 - 0.8e1 * t591 * t261 + 0.8e1 * t591 * t260 * t20;
  t607 = t1 * t150;
  t614 = 0.8e1 * (-t53 + t226 + t55 - t56 + t57 - t566 + 0.16e2) * s * t393 + 0.8e1 * (-t53 + t226 + t56 + t51 - t561 + 0.16e2) * s * t151 - 0.16e2 * t250 * t607 * t13 + 0.16e2 * t250 * t607 * t20;
  t631 = 0.4e1 * t295 * (t53 + t226 + t55 - t56 + t51 - t561 + 0.16e2) * t185 * t18 + 0.4e1 * t295 * (t53 + t226 + t56 - t51 + 0.16e2) * t185 * t6 - 0.8e1 * t591 * t302 + 0.8e1 * t591 * t301 * t20;
  t633 = (-0.128e3 * t10 * t12 * t423 - 0.128e3 * t10 * t12 * t426 + 0.128e3 * t369 * t5) * t40 + (0.16e2 * t369 * t46 * t48 - 0.16e2 * t423 * t60 * t61 - 0.16e2 * t426 * t60 * t61) * t91 + (0.16e2 * t100 * t101 * t423 + 0.16e2 * t100 * t101 * t426 + 0.16e2 * t369 * t94 * t96) * t120 + (-0.32e2 * t123 * t124 * t369 + 0.64e2 * t124 * t42 * t423 + 0.64e2 * t124 * t42 * t426) * t147 + (0.64e2 * t149 * t376 - 0.64e2 * t154 * t423 - 0.64e2 * t154 * t426) * t165 + (-0.64e2 * t167 * t376 - 0.64e2 * t170 * t423 - 0.64e2 * t170 * t426) * t178 + t482 * t209 + t490 * t224 + t500 * t248 + t509 * t277 + t517 * t292 + t526 * t316 + t544 * t326 + t559 * t334 + t578 * t342 + t597 * t350 + t614 * t356 + t631 * t364;
  t636 = t18 * t6;
  t641 = 0.64e2 * t636 - t162 - t163;
  t644 = t151 * t18;
  t651 = -0.64e2 * t150 * t250 + 0.128e3 * t181 * t644 - 0.64e2 * t203 * t250;
  t655 = t11 * t51;
  t659 = -0.16e2 * t33 * t655 + 0.32e2 * t380 * t636 - 0.16e2 * t655;
  t663 = t590 * t150;
  t664 = t663 * t33;
  t669 = 0.32e2 * t252 * t644 - 0.16e2 * t404 * t663 - 0.16e2 * t404 * t664;
  t678 = 0.32e2 * t295 * t644 - 0.16e2 * t413 * t663 - 0.16e2 * t413 * t664;
  t739 = (-0.128e3 * t10 * t12 * t318 - 0.128e3 * t10 * t12 * t321) * t40 + (-0.16e2 * t318 * t60 * t61 - 0.16e2 * t321 * t60 * t61) * t91 + (0.16e2 * t100 * t101 * t318 + 0.16e2 * t100 * t101 * t321) * t120 + (0.64e2 * t124 * t318 * t42 + 0.64e2 * t124 * t321 * t42) * t147 + (-0.64e2 * t154 * t318 - 0.64e2 * t154 * t321) * t165 + (-0.64e2 * t170 * t318 - 0.64e2 * t170 * t321) * t178 - t325 * t209 - t333 * t224 - t341 * t248 - t349 * t277 - t355 * t292 - t363 * t316 + t208 * t326 + t223 * t334 + t247 * t342 + t276 * t350 + t291 * t356 + t315 * t364;
  t742 = t18 + t6;
  t746 = -0.64e2 * t742;
  t751 = -t402;
  t773 = -0.128e3 * t181 * t326 * t376 - 0.32e2 * t252 * t350 * t376 - 0.32e2 * t295 * t364 * t376 - 0.32e2 * t342 * t369 * t380 - 0.32e2 * t356 * t369 * t380 + t120 * t746 - t209 * t397 + t248 * t751 - t277 * t410 + t292 * t751 - t316 * t417 + 0.128e3 * t40 * t742 + t746 * t91;
  t776 = 0.2e1 * beta_y;
  t787 = beta_y * t4;
  t795 = 0.32e2 * t42;
  t801 = 0.5e1 * t54;
  t802 = 0.40e2 * t42;
  t807 = t59 * t4;
  t848 = (t28 - s + 0.8e1) * s;
  t854 = t153 * s;
  t866 = t169 * s;
  t891 = (0.32e2 * (-t28 + t30 - t42 + t776 - 0.4e1) * t2 * t4 * t18 + 0.32e2 * (-t28 + t30 + t42 - t776 - 0.4e1) * t2 * t4 * t6 + 0.64e2 * t787 * t14 - 0.64e2 * t787 * t12 * t20) * t40 + (0.4e1 * (t75 + t76 - t77 - t54 + t795 - t57 - t81 - t561 + 0.16e2) * t4 * t48 * t18 + 0.4e1 * (t75 + t53 - t77 - t801 + t802 + t57 - t81 - t163 + 0.48e2) * t4 * t49 + 0.8e1 * t807 * t48 * t13 - 0.8e1 * t807 * t68) * t91 + (-0.4e1 * (t75 - t53 + t77 - t801 + t802 - t57 - t81 + t163 - 0.48e2) * t2 * t96 * t18 - 0.4e1 * (t75 - t76 + t77 - t54 + t795 + t57 - t81 + t561 - 0.16e2) * t2 * t97 - 0.8e1 * t100 * t96 * t13 + 0.8e1 * t100 * t108) * t120 + (0.16e2 * (t28 - t42 - t528 + 0.4e1) * s * t124 * t18 + 0.16e2 * (t28 + t42 - t528 + 0.4e1) * s * t125 - 0.32e2 * t54 * t128 + 0.32e2 * t54 * t124 * t20) * t147 + (0.16e2 * t151 * t251 * t51 + 0.32e2 * t213 * t854 + 0.16e2 * t393 * t848 - 0.32e2 * t556 * t854) * t165 + (0.16e2 * t294 * t393 * t51 + 0.16e2 * t151 * t848 + 0.32e2 * t213 * t866 - 0.32e2 * t556 * t866) * t178 - t544 * t209 - t559 * t224 - t578 * t248 - t597 * t277 - t614 * t292 - t631 * t316 + t482 * t326 + t490 * t334 + t500 * t342 + t509 * t350 + t517 * t356 + t526 * t364;
  return(At_fH_re * t366 * PREF_V_CF + At_fA_re * t419 * PREF_V_CF + Bt_fH_re * t633 * PREF_V_CF + Bt_fA_re * ((-0.128e3 * t636 + t221 + t222) * t40 + t641 * t91 + t641 * t120 + t651 * t209 + t659 * t248 + t669 * t277 + t659 * t292 + t678 * t316) * PREF_V_CF + At_fH_im * t739 * PREF_V_CF + At_fA_im * t773 * PREF_V_CF + Bt_fH_im * t891 * PREF_V_CF + Bt_fA_im * (t326 * t651 + t342 * t659 + t350 * t669 + t356 * t659 + t364 * t678) * PREF_V_CF);
}
