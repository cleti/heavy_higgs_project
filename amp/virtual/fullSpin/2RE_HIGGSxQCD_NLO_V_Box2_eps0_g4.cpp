
#include "AMP_HEADER.h"

double Eval_V_B2 (AMP_ARGS)
{

AMP_DEFINITIONS

  c_double& cg  = I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_0;
  c_double& cg1 = I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_1;
  c_double& cg3 = I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_0;
  c_double& cg5 = I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_1;


  double t1;
  double t10;
  double t1009;
  double t1014;
  double t102;
  double t1026;
  double t103;
  double t1037;
  double t105;
  double t1055;
  double t1056;
  double t106;
  double t107;
  double t1070;
  double t109;
  double t11;
  double t111;
  double t112;
  double t114;
  double t116;
  double t117;
  double t118;
  double t119;
  double t12;
  double t120;
  double t121;
  double t122;
  double t123;
  double t124;
  double t125;
  double t126;
  double t13;
  double t133;
  double t134;
  double t137;
  double t139;
  double t14;
  double t140;
  double t141;
  double t142;
  double t143;
  double t144;
  double t145;
  double t146;
  double t147;
  double t154;
  double t157;
  double t159;
  double t160;
  double t162;
  double t163;
  double t164;
  double t171;
  double t175;
  double t176;
  double t18;
  double t183;
  double t186;
  double t187;
  double t188;
  double t19;
  double t191;
  double t199;
  double t2;
  double t20;
  double t200;
  double t201;
  double t21;
  double t211;
  double t213;
  double t214;
  double t215;
  double t218;
  double t219;
  double t22;
  double t229;
  double t23;
  double t230;
  double t235;
  double t236;
  double t238;
  double t239;
  double t24;
  double t240;
  double t248;
  double t25;
  double t253;
  double t255;
  double t256;
  double t257;
  double t258;
  double t26;
  double t261;
  double t262;
  double t269;
  double t27;
  double t272;
  double t273;
  double t274;
  double t280;
  double t282;
  double t283;
  double t29;
  double t291;
  double t296;
  double t298;
  double t299;
  double t3;
  double t30;
  double t302;
  double t310;
  double t311;
  double t318;
  double t32;
  double t321;
  double t328;
  double t33;
  double t330;
  double t335;
  double t336;
  double t338;
  double t339;
  double t341;
  double t345;
  double t346;
  double t348;
  double t349;
  double t35;
  double t350;
  double t352;
  double t353;
  double t36;
  double t360;
  double t365;
  double t366;
  double t368;
  double t371;
  double t378;
  double t38;
  double t383;
  double t384;
  double t386;
  double t389;
  double t39;
  double t392;
  double t393;
  double t394;
  double t40;
  double t401;
  double t402;
  double t404;
  double t405;
  double t407;
  double t409;
  double t41;
  double t412;
  double t413;
  double t418;
  double t419;
  double t421;
  double t426;
  double t427;
  double t428;
  double t429;
  double t430;
  double t434;
  double t435;
  double t436;
  double t441;
  double t442;
  double t449;
  double t452;
  double t458;
  double t462;
  double t463;
  double t468;
  double t47;
  double t476;
  double t480;
  double t481;
  double t482;
  double t484;
  double t486;
  double t489;
  double t49;
  double t491;
  double t496;
  double t498;
  double t504;
  double t509;
  double t51;
  double t519;
  double t52;
  double t522;
  double t524;
  double t53;
  double t562;
  double t569;
  double t571;
  double t573;
  double t594;
  double t596;
  double t598;
  double t6;
  double t600;
  double t602;
  double t605;
  double t612;
  double t616;
  double t623;
  double t626;
  double t628;
  double t63;
  double t632;
  double t639;
  double t645;
  double t657;
  double t659;
  double t662;
  double t663;
  double t664;
  double t667;
  double t668;
  double t670;
  double t673;
  double t674;
  double t675;
  double t678;
  double t681;
  double t685;
  double t688;
  double t695;
  double t697;
  double t7;
  double t701;
  double t702;
  double t705;
  double t71;
  double t713;
  double t717;
  double t72;
  double t725;
  double t733;
  double t737;
  double t744;
  double t746;
  double t768;
  double t772;
  double t777;
  double t782;
  double t787;
  double t789;
  double t79;
  double t8;
  double t81;
  double t831;
  double t833;
  double t84;
  double t85;
  double t86;
  double t861;
  double t887;
  double t9;
  double t90;
  double t900;
  double t910;
  double t915;
  double t92;
  double t920;
  double t926;
  double t927;
  double t928;
  double t930;
  double t931;
  double t94;
  double t943;
  double t954;
  double t96;
  double t960;
  double t961;
  double t966;
  double t97;
  double t983;
  double t991;
  double t995;
  t1 = beta * beta;
  t2 = y * y;
  t3 = t1 * t2;
  t6 = beta_y - 0.1e1;
  t7 = 0.1e1 / t6;
  t8 = (0.2e1 * t3 - s + 0.2e1) * t7;
  t9 = beta_y + 0.1e1;
  t10 = 0.1e1 / t9;
  t11 = -s + 0.4e1;
  t12 = 0.1e1 / t11;
  t13 = t10 * t12;
  t14 = sp(s1, k2);
  t18 = 0.1e1 / s;
  t19 = beta_y * t18;
  t20 = t6 * t6;
  t21 = 0.1e1 / t20;
  t22 = t19 * t21;
  t23 = t9 * t9;
  t24 = 0.1e1 / t23;
  t25 = t24 * t12;
  t26 = t1 * t1;
  t27 = s * s;
  t29 = t2 * t2;
  t30 = t26 * t27 * t29;
  t32 = t1 * t27 * t2;
  t33 = 0.4e1 * t32;
  t35 = t1 * s * t2;
  t36 = 0.24e2 * t35;
  t38 = 0.3e1 * t27;
  t39 = 0.8e1 * s;
  t40 = t30 - t33 + t36 - 0.64e2 * t3 + t38 - t39;
  t41 = sp(s1, p1);
  t47 = sp(k1, s2);
  t49 = t21 * t24;
  t51 = t12 * t40;
  t52 = sp(s2, p1);
  t53 = t14 * t52;
  t63 = t26 * s * t29;
  t71 = (t1 * t2 * t26 * t27 * t29 - 0.32e2 * t26 * t29 + 0.24e2 * s - 0.4e1 * t27 - 0.4e1 * t30 + 0.7e1 * t32 - t36 + 0.16e2 * t63 - 0.32e2) * t21;
  t72 = sp(s1, s2);
  t79 = RE(I2_MT2_0_MT2_MU2_0);
  t81 = s * t12;
  t84 = beta_y * s;
  t85 = t12 * t41;
  t86 = t84 * t85;
  t90 = t12 * t14;
  t92 = t84 * t90 * t52;
  t94 = 0.2e1 * s;
  t96 = (t35 - t94 + 0.8e1) * s;
  t97 = t12 * t72;
  t102 = (-0.64e2 * t14 * t81 - 0.32e2 * t86) * t47 + 0.32e2 * t92 - 0.16e2 * t96 * t97 + 0.16e2 * t96 * t12;
  t103 = RE(I2_S12_0_0_MU2_0);
  t105 = 0.2e1 * t84;
  t106 = 0.4e1 * beta_y;
  t107 = 0.3e1 * s;
  t109 = (t35 + t105 + t106 - t107 + 0.4e1) * t7;
  t111 = 0.1e1 / (t84 - s + 0.2e1);
  t112 = t111 * t14;
  t114 = t1 * beta;
  t116 = t2 * y;
  t117 = t114 * t27 * t116;
  t118 = 0.3e1 * t32;
  t119 = 0.16e2 * t35;
  t120 = beta_y * t27;
  t121 = 0.3e1 * t120;
  t122 = 0.24e2 * t84;
  t123 = 0.32e2 * beta_y;
  t124 = t117 - t118 + t119 + t121 - t122 - t27 + t123 + t39 - 0.16e2;
  t125 = t124 * t18;
  t126 = t21 * t111;
  t133 = t111 * t52;
  t134 = t133 * t14;
  t137 = 0.2e1 * t117;
  t139 = t114 * s * t116;
  t140 = 0.12e2 * t139;
  t141 = 0.16e2 * t3;
  t142 = 0.6e1 * t120;
  t143 = 0.28e2 * t84;
  t144 = 0.16e2 * beta_y;
  t145 = 0.16e2 * s;
  t146 = t30 - t137 + t140 + t33 - t36 + t141 - t142 + t143 + t38 - t144 - t145 + 0.16e2;
  t147 = t146 * t21;
  t154 = RE(I2_T11_0_MT2_MU2_0);
  t157 = (t35 - t105 - t106 - t107 + 0.4e1) * t10;
  t159 = 0.1e1 / (t84 + s - 0.2e1);
  t160 = t159 * t14;
  t162 = t117 + t118 - t119 + t121 - t122 + t27 + t123 - t39 + 0.16e2;
  t163 = t162 * t18;
  t164 = t24 * t159;
  t171 = t159 * t52;
  t175 = t30 + t137 - t140 + t33 - t36 + t141 + t142 - t143 + t38 + t144 - t145 + 0.16e2;
  t176 = t175 * t24;
  t183 = RE(I2_T12_0_MT2_MU2_0);
  t186 = (t35 - s + 0.2e1) * s;
  t187 = t159 * t111;
  t188 = t187 * t14;
  t191 = t187 * t41;
  t199 = -s + 0.2e1;
  t200 = t35 - s + 0.4e1;
  t201 = t199 * t200;
  t211 = RE(I1_MT2_MU2_0);
  t213 = t199 * t27;
  t214 = t200 * t200;
  t215 = 0.1e1 / t214;
  t218 = t215 * t199;
  t219 = t218 * t41;
  t229 = 0.1e1 / t200;
  t230 = t229 * t72;
  t235 = (-0.128e3 * t14 * t213 * t215 - 0.128e3 * t120 * t219) * t47 + 0.128e3 * t120 * t215 * t199 * t14 * t52 - 0.64e2 * t213 * t230 + 0.64e2 * t213 * t229;
  t236 = RE(I3_0_0_S12_0_0_0_MU2_1);
  t238 = t14 * s;
  t239 = beta_y - 0.2e1;
  t240 = s * t239;
  t248 = s * (t35 - t105 - s + 0.4e1);
  t253 = RE(I3_0_T11_MT2_0_0_MT2_MU2_0);
  t255 = t35 + t84 - t94 + 0.4e1;
  t256 = t20 * t255;
  t257 = t27 * t215;
  t258 = t257 * t14;
  t261 = t20 * (t84 - s + 0.4e1);
  t262 = t257 * t41;
  t269 = t215 * t52 * t14;
  t272 = t27 * t20;
  t273 = t199 * t229;
  t274 = t273 * t72;
  t280 = RE(I3_0_T11_MT2_0_0_MT2_MU2_1);
  t282 = beta_y + 0.2e1;
  t283 = s * t282;
  t291 = s * (t35 + t105 - s + 0.4e1);
  t296 = RE(I3_0_T12_MT2_0_0_MT2_MU2_0);
  t298 = t35 - t84 - t94 + 0.4e1;
  t299 = t23 * t298;
  t302 = t23 * (t84 + s - 0.4e1);
  t310 = t23 * t199;
  t311 = t27 * t229;
  t318 = RE(I3_0_T12_MT2_0_0_MT2_MU2_1);
  t321 = (0.4e1 - t107) * s;
  t328 = 0.4e1 * t35;
  t330 = (t328 - t38 + t145 - 0.16e2) * s;
  t335 = (0.16e2 * t321 * t90 - 0.64e2 * t86) * t47 + 0.64e2 * t92 - 0.8e1 * t330 * t97 + 0.8e1 * t330 * t12;
  t336 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t338 = t47 * t14;
  t339 = t338 * t27;
  t341 = t11 * t27;
  t345 = -0.8e1 * t341 * t72 - 0.16e2 * t339 + 0.8e1 * t341;
  t346 = RE(cg);
  t348 = t27 * s;
  t349 = t348 * t215;
  t350 = t349 * t14;
  t352 = t349 * t41;
  t353 = t261 * t352;
  t360 = t348 * t20;
  t365 = (-0.16e2 * t256 * t350 + 0.16e2 * t353) * t47 - 0.16e2 * t261 * t348 * t269 - 0.16e2 * t360 * t274 + 0.16e2 * t360 * t273;
  t366 = RE(cg1);
  t368 = RE(cg3);
  t371 = t302 * t352;
  t378 = t348 * t23;
  t383 = (-0.16e2 * t299 * t350 + 0.16e2 * t371) * t47 - 0.16e2 * t302 * t348 * t269 - 0.16e2 * t378 * t274 + 0.16e2 * t378 * t273;
  t384 = RE(cg5);
  t386 = EPS_(k1, k2, p1, s1);
  t389 = EPS_(k1, k2, p1, s2);
  t392 = t12 * t386 * t84 + t12 * t389 * t84;
  t393 = 0.32e2 * t392;
  t394 = IM(I2_S12_0_0_MU2_0);
  t401 = 0.128e3 * t120 * t218 * t386 + 0.128e3 * t120 * t218 * t389;
  t402 = IM(I3_0_0_S12_0_0_0_MU2_1);
  t404 = 0.64e2 * t392;
  t405 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t407 = t349 * t386;
  t409 = t349 * t389;
  t412 = -0.16e2 * t261 * t407 - 0.16e2 * t261 * t409;
  t413 = IM(cg1);
  t418 = -0.16e2 * t302 * t407 - 0.16e2 * t302 * t409;
  t419 = IM(cg5);
  t421 = ((0.32e2 * t22 * t25 * t40 * t41 + 0.128e3 * t13 * t14 * t8) * t47 - 0.32e2 * t19 * t49 * t51 * t53 + 0.16e2 * t71 * t25 * t72 - 0.16e2 * t71 * t25) * t79 + t102 * t103 + ((0.16e2 * t125 * t126 * t41 - 0.16e2 * t109 * t112) * t47 - 0.16e2 * t125 * t21 * t134 + 0.8e1 * t147 * t111 * t72 - 0.8e1 * t147 * t111) * t154 + ((-0.16e2 * t163 * t164 * t41 - 0.16e2 * t157 * t160) * t47 + 0.16e2 * t163 * t24 * t171 * t14 - 0.8e1 * t176 * t159 * t72 + 0.8e1 * t176 * t159) * t183 + ((-0.32e2 * t186 * t188 + 0.64e2 * t191 * t84) * t47 - 0.64e2 * t84 * t159 * t134 - 0.16e2 * t201 * s * t187 * t72 + 0.16e2 * t201 * s * t159 * t111) * t211 + t235 * t236 + ((0.16e2 * t240 * t41 + 0.16e2 * t238) * t47 - 0.16e2 * t240 * t53 + 0.8e1 * t248 * t72 - 0.8e1 * t248) * t253 + ((0.32e2 * t256 * t258 - 0.32e2 * t261 * t262) * t47 + 0.32e2 * t261 * t27 * t269 + 0.32e2 * t272 * t274 - 0.32e2 * t272 * t273) * t280 + ((0.16e2 * t283 * t41 + 0.16e2 * t238) * t47 - 0.16e2 * t283 * t53 + 0.8e1 * t291 * t72 - 0.8e1 * t291) * t296 + ((0.32e2 * t258 * t299 - 0.32e2 * t262 * t302) * t47 + 0.32e2 * t302 * t27 * t269 + 0.32e2 * t310 * t311 * t72 - 0.32e2 * t310 * t311) * t318 + t335 * t336 + t345 * t346 + t365 * t366 + t345 * t368 + t383 * t384 + t393 * t394 + t401 * t402 + t404 * t405 + t412 * t413 + t418 * t419;
  t426 = t63 - 0.3e1 * t35 - 0.4e1 * t3 + t94 - 0.4e1;
  t427 = t426 * t18;
  t428 = t427 * t21;
  t429 = EPS_(k1, k2, s1, s2);
  t430 = t24 * t429;
  t434 = 0.4e1 * t84;
  t435 = t35 - t434 + t107 - 0.4e1;
  t436 = t435 * t18;
  t441 = t35 + t434 + t107 - 0.4e1;
  t442 = t441 * t18;
  t449 = s * t429;
  t452 = t229 * t429;
  t458 = t27 * t23;
  t462 = t27 * t429;
  t463 = t462 * t346;
  t468 = t462 * t368;
  t476 = -0.256e3 * t14 * t311 - 0.256e3 * t311 * t47;
  t480 = -t14 * t27 - t27 * t47;
  t481 = 0.32e2 * t480;
  t482 = IM(cg);
  t484 = t229 * t47;
  t486 = t229 * t14;
  t489 = -0.64e2 * t360 * t484 - 0.64e2 * t360 * t486;
  t491 = IM(cg3);
  t496 = -0.64e2 * t378 * t484 - 0.64e2 * t378 * t486;
  t498 = 0.128e3 * t428 * t430 * t79 - 0.64e2 * t436 * t21 * t429 * t154 - 0.64e2 * t442 * t430 * t183 + 0.256e3 * t311 * t429 * t236 - 0.64e2 * t449 * t253 - 0.128e3 * t272 * t452 * t280 - 0.64e2 * t449 * t296 - 0.128e3 * t458 * t452 * t318 + 0.32e2 * t463 + 0.64e2 * t360 * t452 * t366 + 0.32e2 * t468 + 0.64e2 * t378 * t452 * t384 + t476 * t402 + t481 * t482 + t489 * t413 + t481 * t491 + t496 * t419;
  t504 = EPS_(k1, p1, s1, s2);
  t509 = EPS_(k2, p1, s1, s2);
  t519 = t84 * t12 * t504;
  t522 = t84 * t12 * t509;
  t524 = -0.64e2 * t429 * t81 - 0.32e2 * t519 - 0.32e2 * t522;
  t562 = -0.128e3 * t120 * t218 * t504 - 0.128e3 * t120 * t218 * t509 - 0.128e3 * t213 * t215 * t429;
  t569 = t257 * t429;
  t571 = t257 * t504;
  t573 = t257 * t509;
  t594 = 0.16e2 * t12 * t321 * t429 - 0.64e2 * t519 - 0.64e2 * t522;
  t596 = (0.32e2 * t22 * t25 * t40 * t504 + 0.32e2 * t22 * t25 * t40 * t509 + 0.128e3 * t13 * t429 * t8) * t79 + t524 * t103 + (-0.16e2 * t109 * t111 * t429 + 0.16e2 * t125 * t126 * t504 + 0.16e2 * t125 * t126 * t509) * t154 + (-0.16e2 * t157 * t159 * t429 - 0.16e2 * t163 * t164 * t504 - 0.16e2 * t163 * t164 * t509) * t183 + (-0.32e2 * t186 * t187 * t429 + 0.64e2 * t187 * t504 * t84 + 0.64e2 * t187 * t509 * t84) * t211 + t562 * t236 + (0.16e2 * t240 * t504 + 0.16e2 * t240 * t509 + 0.16e2 * t449) * t253 + (0.32e2 * t256 * t569 - 0.32e2 * t261 * t571 - 0.32e2 * t261 * t573) * t280 + (0.16e2 * t283 * t504 + 0.16e2 * t283 * t509 + 0.16e2 * t449) * t296 + (0.32e2 * t299 * t569 - 0.32e2 * t302 * t571 - 0.32e2 * t302 * t573) * t318 + t594 * t336;
  t598 = t349 * t429;
  t600 = t349 * t504;
  t602 = t349 * t509;
  t605 = -0.16e2 * t256 * t598 + 0.16e2 * t261 * t600 + 0.16e2 * t261 * t602;
  t612 = -0.16e2 * t299 * t598 + 0.16e2 * t302 * t600 + 0.16e2 * t302 * t602;
  t616 = t12 * t47;
  t623 = t120 * t85;
  t626 = t120 * t12 * t52;
  t628 = -0.8e1 * (t35 + t84 + 0.8e1) * s * t616 - 0.8e1 * (t35 - t84 + 0.8e1) * s * t90 - 0.16e2 * t623 + 0.16e2 * t626;
  t632 = t257 * t47;
  t639 = beta_y * t348;
  t645 = -0.32e2 * t199 * (t35 + t84 + 0.4e1) * t632 - 0.32e2 * t199 * (t35 - t84 + 0.4e1) * t258 - 0.64e2 * t639 * t219 + 0.64e2 * t639 * t218 * t52;
  t657 = -0.16e2 * (t35 + t84 + t107 - 0.4e1) * s * t616 - 0.16e2 * (t35 - t84 + t107 - 0.4e1) * s * t90 - 0.32e2 * t623 + 0.32e2 * t626;
  t659 = 0.16e2 * t480;
  t662 = -t32 + t328 + t27 - 0.12e2 * s + 0.16e2;
  t663 = t20 * t662;
  t664 = t349 * t47;
  t667 = 0.2e1 * t120;
  t668 = 0.4e1 * s;
  t670 = t20 * (-t32 + t328 + t667 - t27 - t668 + 0.16e2);
  t673 = t27 * t27;
  t674 = t673 * t215;
  t675 = t674 * t41;
  t678 = t674 * t52;
  t681 = 0.8e1 * t261 * t675 - 0.8e1 * t261 * t678 - 0.4e1 * t350 * t670 - 0.4e1 * t663 * t664;
  t685 = t23 * (-t32 + t328 - t667 - t27 - t668 + 0.16e2);
  t688 = t23 * t662;
  t695 = 0.8e1 * t302 * t675 - 0.8e1 * t302 * t678 - 0.4e1 * t350 * t688 - 0.4e1 * t664 * t685;
  t697 = t366 * t605 + t384 * t612 + t394 * t628 + t402 * t645 + t405 * t657 + t413 * t681 + t419 * t695 + t482 * t659 + t491 * t659 - 0.16e2 * t463 - 0.16e2 * t468;
  t701 = t24 * t14;
  t702 = t701 * t47;
  t705 = t426 * t21;
  t713 = t21 * t14;
  t717 = t435 * t21;
  t725 = t441 * t24;
  t733 = t348 * t229;
  t737 = -0.256e3 * t311 * t338 + 0.128e3 * t72 * t733 + 0.128e3 * t733;
  t744 = 0.64e2 * s * t338 - 0.32e2 * t27 * t72 - 0.32e2 * t27;
  t746 = t486 * t47;
  t768 = 0.16e2 * t348 * t72 - 0.32e2 * t339 + 0.16e2 * t348;
  t772 = t673 * t20;
  t777 = 0.32e2 * t229 * t772 + 0.32e2 * t230 * t772 - 0.64e2 * t360 * t746;
  t782 = t673 * t23;
  t787 = 0.32e2 * t229 * t782 + 0.32e2 * t230 * t782 - 0.64e2 * t378 * t746;
  t789 = (0.64e2 * t24 * t705 * t72 + 0.64e2 * t24 * t705 - 0.128e3 * t428 * t702) * t79 + (0.64e2 * t436 * t47 * t713 - 0.32e2 * t717 * t72 - 0.32e2 * t717) * t154 + (0.64e2 * t442 * t702 - 0.32e2 * t72 * t725 - 0.32e2 * t725) * t183 + t737 * t236 + t744 * t253 + (-0.64e2 * t229 * t360 - 0.64e2 * t230 * t360 + 0.128e3 * t272 * t746) * t280 + t744 * t296 + (-0.64e2 * t229 * t378 - 0.64e2 * t230 * t378 + 0.128e3 * t458 * t746) * t318 + t768 * t346 + t777 * t366 + t768 * t368 + t787 * t384;
  t831 = t257 * t386;
  t833 = t257 * t389;
  t861 = (0.32e2 * t22 * t25 * t386 * t40 + 0.32e2 * t22 * t25 * t389 * t40) * t79 - t393 * t103 + (0.16e2 * t125 * t126 * t386 + 0.16e2 * t125 * t126 * t389) * t154 + (-0.16e2 * t163 * t164 * t386 - 0.16e2 * t163 * t164 * t389) * t183 + (0.64e2 * t187 * t386 * t84 + 0.64e2 * t187 * t389 * t84) * t211 - t401 * t236 + (0.16e2 * t240 * t386 + 0.16e2 * t240 * t389) * t253 + (-0.32e2 * t261 * t831 - 0.32e2 * t261 * t833) * t280 + (0.16e2 * t283 * t386 + 0.16e2 * t283 * t389) * t296 + (-0.32e2 * t302 * t831 - 0.32e2 * t302 * t833) * t318 - t404 * t336 - t412 * t366 - t418 * t384 + t102 * t394 + t235 * t402 + t335 * t405 + t345 * t482 + t365 * t413 + t345 * t491 + t383 * t419;
  t887 = -0.64e2 * s * t47 - 0.64e2 * t238;
  t900 = -t481;
  t910 = t462 * t482;
  t915 = t462 * t491;
  t920 = (0.128e3 * t14 * t427 * t49 + 0.128e3 * t427 * t47 * t49) * t79 + (-0.64e2 * t21 * t436 * t47 - 0.64e2 * t436 * t713) * t154 + (-0.64e2 * t24 * t442 * t47 - 0.64e2 * t442 * t701) * t183 - t476 * t236 + t887 * t253 + (-0.128e3 * t272 * t484 - 0.128e3 * t272 * t486) * t280 + t887 * t296 + (-0.128e3 * t458 * t484 - 0.128e3 * t458 * t486) * t318 + t900 * t346 - t489 * t366 + t900 * t368 - t496 * t384 + 0.256e3 * t311 * t429 * t402 + 0.32e2 * t910 + 0.64e2 * t360 * t452 * t413 + 0.32e2 * t915 + 0.64e2 * t378 * t452 * t419;
  t926 = t26 * beta * t27 * t29 * y;
  t927 = 0.4e1 * t117;
  t928 = 0.24e2 * t139;
  t930 = 0.32e2 * t114 * t116;
  t931 = 0.32e2 * t3;
  t943 = beta_y * t21 * t24;
  t954 = 0.12e2 * t35;
  t960 = 0.32e2 * t84;
  t961 = 0.20e2 * s;
  t966 = t124 * t21;
  t983 = t162 * t24;
  t991 = t298 * s;
  t995 = t255 * s;
  t1009 = 0.3e1 * t84;
  t1014 = t27 * t239;
  t1026 = t349 * t52;
  t1037 = t27 * t282;
  t1055 = (-0.8e1 * (t926 - t927 + t928 - t930 - t931 + t121 - t122 + t123 + t145 - 0.32e2) * t21 * t13 * t47 - 0.8e1 * (t926 - t927 + t928 - t930 + t931 + t121 - t122 + t123 - t145 + 0.32e2) * t7 * t25 * t14 - 0.16e2 * t943 * t51 * t41 + 0.16e2 * t943 * t51 * t52) * t79 - t628 * t103 + (-0.4e1 * (t30 - t137 + t140 - t954 + t141 + t667 + t434 - t27 + t144 - t668) * t21 * t111 * t47 - 0.4e1 * (t117 - t118 + t954 + t121 - t960 - t27 + t144 + t961 - 0.32e2) * t7 * t112 - 0.8e1 * t966 * t111 * t41 + 0.8e1 * t966 * t133) * t154 + (0.4e1 * (t117 + t118 - t954 + t121 - t960 + t27 + t144 - t961 + 0.32e2) * t10 * t159 * t47 + 0.4e1 * (t30 + t137 - t140 - t954 + t141 - t667 - t434 - t27 - t144 - t668) * t24 * t160 + 0.8e1 * t983 * t159 * t41 - 0.8e1 * t983 * t171) * t183 + (0.32e2 * t120 * t187 * t52 + 0.16e2 * t187 * t47 * t991 - 0.32e2 * t120 * t191 + 0.16e2 * t188 * t995) * t211 - t645 * t236 + (-0.4e1 * t991 * t47 - 0.4e1 * s * (t35 - t1009 + t94 + 0.4e1) * t14 - 0.8e1 * t1014 * t41 + 0.8e1 * t1014 * t52) * t253 + (-0.16e2 * t1026 * t261 - 0.8e1 * t258 * t670 - 0.8e1 * t632 * t663 + 0.16e2 * t353) * t280 + (-0.4e1 * s * (t35 + t1009 + t94 + 0.4e1) * t47 - 0.4e1 * t995 * t14 - 0.8e1 * t1037 * t41 + 0.8e1 * t1037 * t52) * t296 + (-0.16e2 * t1026 * t302 - 0.8e1 * t258 * t688 - 0.8e1 * t632 * t685 + 0.16e2 * t371) * t318 - t657 * t336;
  t1056 = -t659;
  t1070 = t1056 * t346 + t1056 * t368 - t366 * t681 - t384 * t695 + t394 * t524 + t402 * t562 + t405 * t594 + t413 * t605 + t419 * t612 - 0.16e2 * t910 - 0.16e2 * t915;
  return(At_fH_re * t421 * PREF_V_CA + At_fA_re * t498 * PREF_V_CA + Bt_fH_re * (t596 + t697) * PREF_V_CA + Bt_fA_re * t789 * PREF_V_CA + At_fH_im * t861 * PREF_V_CA + At_fA_im * t920 * PREF_V_CA + Bt_fH_im * (t1055 + t1070) * PREF_V_CA + Bt_fA_im * (t402 * t737 + t413 * t777 + t419 * t787 + t482 * t768 + t491 * t768) * PREF_V_CA);
}