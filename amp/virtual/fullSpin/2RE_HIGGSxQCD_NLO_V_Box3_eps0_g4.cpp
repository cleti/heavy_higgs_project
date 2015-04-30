
#include "AMP_HEADER.h"

double Eval_V_B3 (AMP_ARGS)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_V(ap);

  c_double& cg  = I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_0;
  c_double& cg1 = I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_1;
  // c_double& cg3 = I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_0;
  // c_double& cg5 = I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_1;


  double t1;
  double t10;
  double t100;
  double t101;
  double t1013;
  double t102;
  double t1026;
  double t103;
  double t104;
  double t111;
  double t114;
  double t116;
  double t117;
  double t119;
  double t120;
  double t121;
  double t128;
  double t13;
  double t132;
  double t133;
  double t140;
  double t143;
  double t144;
  double t145;
  double t148;
  double t15;
  double t157;
  double t158;
  double t16;
  double t168;
  double t17;
  double t170;
  double t171;
  double t172;
  double t174;
  double t175;
  double t176;
  double t179;
  double t18;
  double t180;
  double t181;
  double t182;
  double t183;
  double t184;
  double t189;
  double t19;
  double t190;
  double t193;
  double t195;
  double t2;
  double t200;
  double t202;
  double t203;
  double t205;
  double t206;
  double t207;
  double t208;
  double t21;
  double t211;
  double t212;
  double t213;
  double t22;
  double t221;
  double t222;
  double t224;
  double t225;
  double t229;
  double t23;
  double t233;
  double t236;
  double t239;
  double t24;
  double t240;
  double t248;
  double t253;
  double t255;
  double t256;
  double t259;
  double t26;
  double t260;
  double t261;
  double t269;
  double t27;
  double t276;
  double t278;
  double t279;
  double t28;
  double t280;
  double t281;
  double t284;
  double t285;
  double t293;
  double t298;
  double t3;
  double t300;
  double t303;
  double t304;
  double t312;
  double t317;
  double t32;
  double t321;
  double t324;
  double t325;
  double t326;
  double t335;
  double t337;
  double t338;
  double t34;
  double t340;
  double t342;
  double t345;
  double t346;
  double t347;
  double t348;
  double t35;
  double t356;
  double t358;
  double t359;
  double t364;
  double t365;
  double t367;
  double t370;
  double t374;
  double t375;
  double t382;
  double t383;
  double t385;
  double t388;
  double t389;
  double t39;
  double t390;
  double t391;
  double t392;
  double t393;
  double t397;
  double t398;
  double t403;
  double t404;
  double t408;
  double t409;
  double t41;
  double t410;
  double t414;
  double t418;
  double t419;
  double t42;
  double t426;
  double t43;
  double t430;
  double t435;
  double t436;
  double t446;
  double t453;
  double t455;
  double t46;
  double t460;
  double t463;
  double t48;
  double t49;
  double t5;
  double t50;
  double t500;
  double t503;
  double t508;
  double t537;
  double t56;
  double t562;
  double t572;
  double t576;
  double t58;
  double t583;
  double t589;
  double t59;
  double t593;
  double t598;
  double t6;
  double t60;
  double t601;
  double t603;
  double t609;
  double t611;
  double t614;
  double t615;
  double t618;
  double t62;
  double t625;
  double t629;
  double t637;
  double t64;
  double t643;
  double t646;
  double t65;
  double t652;
  double t663;
  double t67;
  double t69;
  double t693;
  double t697;
  double t7;
  double t70;
  double t701;
  double t707;
  double t709;
  double t71;
  double t72;
  double t73;
  double t74;
  double t740;
  double t742;
  double t75;
  double t76;
  double t77;
  double t782;
  double t785;
  double t8;
  double t810;
  double t812;
  double t84;
  double t847;
  double t85;
  double t851;
  double t861;
  double t869;
  double t870;
  double t871;
  double t872;
  double t877;
  double t878;
  double t88;
  double t883;
  double t896;
  double t9;
  double t90;
  double t900;
  double t91;
  double t925;
  double t926;
  double t927;
  double t928;
  double t929;
  double t93;
  double t934;
  double t935;
  double t936;
  double t937;
  double t938;
  double t94;
  double t943;
  double t95;
  double t96;
  double t969;
  double t97;
  double t98;
  double t980;
  double t988;
  double t99;
  double t997;
  t1 = beta * beta;
  t2 = y * y;
  t3 = t1 * t2;
  t5 = beta_y - 0.1e1;
  t6 = 0.1e1 / t5;
  t7 = (t3 + 0.3e1) * t6;
  t8 = beta_y + 0.1e1;
  t9 = 0.1e1 / t8;
  t10 = sp(s1, k2);
  t13 = 0.1e1 / s;
  t15 = t5 * t5;
  t16 = 0.1e1 / t15;
  t17 = beta_y * t13 * t16;
  t18 = t8 * t8;
  t19 = 0.1e1 / t18;
  t21 = t1 * s * t2;
  t22 = 0.3e1 * t21;
  t23 = 0.12e2 * t3;
  t24 = 0.3e1 * s;
  t26 = t19 * (-t22 + t23 + t24 - 0.4e1);
  t27 = sp(s1, p1);
  t28 = t26 * t27;
  t32 = sp(k1, s2);
  t34 = sp(s2, p1);
  t35 = t10 * t34;
  t39 = t1 * t1;
  t41 = t2 * t2;
  t42 = t39 * s * t41;
  t43 = 0.2e1 * t42;
  t46 = 0.5e1 * t21;
  t48 = (0.8e1 * t39 * t41 - t23 - t24 - t43 + t46 + 0.12e2) * t16;
  t49 = sp(s1, s2);
  t50 = t19 * t49;
  t56 = RE(I2_MT2_0_MT2_MU2_0);
  t58 = beta_y * s;
  t59 = 0.2e1 * t58;
  t60 = 0.2e1 * s;
  t62 = (t59 + beta_y - t60 + 0.3e1) * t6;
  t64 = 0.1e1 / (t58 - s + 0.2e1);
  t65 = t64 * t10;
  t67 = s * s;
  t69 = t1 * t67 * t2;
  t70 = beta_y * t67;
  t71 = 0.2e1 * t70;
  t72 = 0.10e2 * t58;
  t73 = 0.12e2 * beta_y;
  t74 = 0.5e1 * s;
  t75 = -t69 + t46 + t71 - t72 - t67 + t73 + t74 - 0.8e1;
  t76 = t75 * t13;
  t77 = t16 * t64;
  t84 = t64 * t34;
  t85 = t84 * t10;
  t88 = t1 * beta;
  t90 = t2 * y;
  t91 = t88 * t67 * t90;
  t93 = t88 * s * t90;
  t94 = 0.5e1 * t93;
  t95 = 0.4e1 * t69;
  t96 = 0.17e2 * t21;
  t97 = 0.8e1 * t3;
  t98 = 0.5e1 * t70;
  t99 = 0.23e2 * t58;
  t100 = 0.2e1 * t67;
  t101 = 0.16e2 * beta_y;
  t102 = 0.11e2 * s;
  t103 = -t91 + t94 + t95 - t96 + t97 - t98 + t99 + t100 - t101 - t102 + 0.12e2;
  t104 = t103 * t16;
  t111 = RE(I2_T11_0_MT2_MU2_0);
  t114 = (t59 + beta_y + t60 - 0.3e1) * t9;
  t116 = 0.1e1 / (t58 + s - 0.2e1);
  t117 = t116 * t10;
  t119 = -t69 + t46 - t71 + t72 - t67 - t73 + t74 - 0.8e1;
  t120 = t119 * t13;
  t121 = t19 * t116;
  t128 = t116 * t34;
  t132 = -t91 + t94 - t95 + t96 - t97 - t98 + t99 - t100 - t101 + t102 - 0.12e2;
  t133 = t132 * t19;
  t140 = RE(I2_T12_0_MT2_MU2_0);
  t143 = (t21 - s + 0.2e1) * s;
  t144 = t116 * t64;
  t145 = t144 * t10;
  t148 = t144 * t27;
  t157 = t21 - s + 0.4e1;
  t158 = (-s + 0.2e1) * t157;
  t168 = RE(I1_MT2_MU2_0);
  t170 = 0.3e1 * t93;
  t171 = 0.3e1 * t58;
  t172 = 0.8e1 * beta_y;
  t174 = (t170 - t21 - t171 + t172 + s) * s;
  t175 = 0.1e1 / t157;
  t176 = t175 * t10;
  t179 = 0.2e1 * t21;
  t180 = 0.2e1 * t3;
  t181 = 0.6e1 * beta_y;
  t182 = t93 - t179 - t180 - t58 + t181 + t60 - 0.8e1;
  t183 = t182 * s;
  t184 = t175 * t27;
  t189 = t175 * t34;
  t190 = t189 * t10;
  t193 = 0.7e1 * t58;
  t195 = (t179 - t193 + t172 + s) * s;
  t200 = RE(I3_0_T11_MT2_0_0_MT2_MU2_0);
  t202 = t67 * t8;
  t203 = t202 * t15;
  t205 = t157 * t157;
  t206 = 0.1e1 / t205;
  t207 = (t22 - t24 + 0.8e1) * t206;
  t208 = t207 * t10;
  t211 = t70 * t15;
  t212 = t206 * t8;
  t213 = t212 * t27;
  t221 = t8 * t15;
  t222 = 0.8e1 - t24;
  t224 = t67 * t175;
  t225 = t224 * t49;
  t229 = t222 * t67 * t175;
  t233 = RE(I3_0_T11_MT2_0_0_MT2_MU2_1);
  t236 = (t170 + t21 - t171 + t172 - s) * s;
  t239 = t93 + t179 + t180 - t58 + t181 - t60 + 0.8e1;
  t240 = t239 * s;
  t248 = (t179 + t193 - t172 + s) * s;
  t253 = RE(I3_0_T12_MT2_0_0_MT2_MU2_0);
  t255 = t67 * t5;
  t256 = t255 * t18;
  t259 = t70 * t18;
  t260 = t206 * t5;
  t261 = t260 * t27;
  t269 = t5 * t18;
  t276 = RE(I3_0_T12_MT2_0_0_MT2_MU2_1);
  t278 = t22 - t24 + 0.16e2;
  t279 = t5 * t278;
  t280 = s * t175;
  t281 = t280 * t10;
  t284 = t21 + t58 - t60 + 0.8e1;
  t285 = t284 * t175;
  t293 = (-t171 + t172 + t24 - 0.16e2) * s;
  t298 = RE(I3_MT2_0_T11_0_MT2_MT2_MU2_0);
  t300 = t8 * t278;
  t303 = t21 - t58 - t60 + 0.8e1;
  t304 = t303 * t175;
  t312 = (-t171 + t172 - t24 + 0.16e2) * s;
  t317 = RE(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
  t321 = (0.3e1 * t42 + t179 + t97 - t74 + 0.24e2) * t67;
  t324 = t70 * t175;
  t325 = t5 * t8;
  t326 = t325 * t27;
  t335 = (-t22 + t97 - t74 + 0.24e2) * t67;
  t337 = (-0.2e1 * t176 * t321 + 0.8e1 * t324 * t326) * t32 - 0.8e1 * t324 * t325 * t35 - t335 * t49 + t335;
  t338 = RE(cg);
  t340 = t67 * s;
  t342 = t340 * t15 * t18;
  t345 = beta_y * t340;
  t346 = t345 * t15;
  t347 = t18 * t206;
  t348 = t347 * t27;
  t356 = t15 * t18;
  t358 = t340 * t175;
  t359 = t358 * t49;
  t364 = (-0.2e1 * t208 * t342 + 0.8e1 * t346 * t348) * t32 - 0.8e1 * t346 * t347 * t35 - t356 * t222 * t359 + t356 * t222 * t340 * t175;
  t365 = RE(cg1);
  t367 = EPS_(k1, k2, p1, s1);
  t370 = EPS_(k1, k2, p1, s2);
  t374 = -0.8e1 * t324 * t325 * t367 - 0.8e1 * t324 * t325 * t370;
  t375 = IM(cg);
  t382 = -0.8e1 * t346 * t347 * t367 - 0.8e1 * t346 * t347 * t370;
  t383 = IM(cg1);
  t385 = ((-0.32e2 * t10 * t7 * t9 + 0.32e2 * t17 * t28) * t32 - 0.32e2 * t17 * t26 * t35 + 0.16e2 * t48 * t50 - 0.16e2 * t48 * t19) * t56 + ((-0.32e2 * t27 * t76 * t77 + 0.32e2 * t62 * t65) * t32 + 0.32e2 * t76 * t16 * t85 - 0.16e2 * t104 * t64 * t49 + 0.16e2 * t104 * t64) * t111 + ((-0.32e2 * t120 * t121 * t27 - 0.32e2 * t114 * t117) * t32 + 0.32e2 * t120 * t19 * t128 * t10 - 0.16e2 * t133 * t116 * t49 + 0.16e2 * t133 * t116) * t140 + ((0.32e2 * t143 * t145 - 0.64e2 * t148 * t58) * t32 + 0.64e2 * t58 * t116 * t85 + 0.16e2 * t158 * s * t144 * t49 - 0.16e2 * t158 * s * t116 * t64) * t168 + ((-0.4e1 * t174 * t176 - 0.8e1 * t183 * t184) * t32 + 0.8e1 * t183 * t190 - 0.2e1 * t195 * t49 + 0.2e1 * t195) * t200 + ((-0.4e1 * t203 * t208 + 0.16e2 * t211 * t213) * t32 - 0.16e2 * t211 * t212 * t35 - 0.2e1 * t221 * t222 * t225 + 0.2e1 * t221 * t229) * t233 + ((0.4e1 * t176 * t236 - 0.8e1 * t184 * t240) * t32 + 0.8e1 * t240 * t190 - 0.2e1 * t248 * t49 + 0.2e1 * t248) * t253 + ((0.4e1 * t208 * t256 - 0.16e2 * t259 * t261) * t32 + 0.16e2 * t259 * t260 * t35 + 0.2e1 * t269 * t222 * t225 - 0.2e1 * t269 * t229) * t276 + ((0.16e2 * t27 * t285 - 0.4e1 * t279 * t281) * t32 - 0.16e2 * t285 * t35 - 0.2e1 * t293 * t49 + 0.2e1 * t293) * t298 + ((-0.16e2 * t27 * t304 + 0.4e1 * t281 * t300) * t32 + 0.16e2 * t304 * t35 + 0.2e1 * t312 * t49 - 0.2e1 * t312) * t317 + t337 * t338 + t364 * t365 + t374 * t375 + t382 * t383;
  t388 = 0.4e1 * t3;
  t389 = t43 - t46 - t388 + t24 - 0.4e1;
  t390 = t389 * t13;
  t391 = t390 * t16;
  t392 = EPS_(k1, k2, s1, s2);
  t393 = t19 * t392;
  t397 = t21 - t171 + t60 - 0.2e1;
  t398 = t397 * t13;
  t403 = t21 + t171 + t60 - 0.2e1;
  t404 = t403 * t13;
  t408 = 0.3e1 * beta_y;
  t409 = t408 + 0.1e1;
  t410 = t409 * s;
  t414 = t175 * t392;
  t418 = t408 - 0.1e1;
  t419 = t418 * s;
  t426 = t5 * s;
  t430 = t8 * s;
  t435 = 0.3e1 * t3 + 0.5e1;
  t436 = t435 * t67;
  t446 = -0.4e1 * t10 * t436 - 0.4e1 * t32 * t436;
  t453 = -0.12e2 * t10 * t356 * t358 - 0.12e2 * t32 * t356 * t358;
  t455 = -0.64e2 * t391 * t393 * t56 + 0.64e2 * t398 * t16 * t392 * t111 + 0.64e2 * t404 * t393 * t140 + 0.8e1 * t410 * t392 * t200 + 0.24e2 * t203 * t414 * t233 - 0.8e1 * t419 * t392 * t253 - 0.24e2 * t256 * t414 * t276 + 0.24e2 * t426 * t392 * t298 - 0.24e2 * t430 * t392 * t317 + 0.4e1 * t436 * t392 * t338 + 0.12e2 * t342 * t414 * t365 + t446 * t375 + t453 * t383;
  t460 = EPS_(k1, p1, s1, s2);
  t463 = EPS_(k2, p1, s1, s2);
  t500 = t175 * t460;
  t503 = t175 * t463;
  t508 = t207 * t392;
  t537 = t280 * t392;
  t562 = 0.8e1 * t324 * t325 * t460 + 0.8e1 * t324 * t325 * t463 - 0.2e1 * t321 * t414;
  t572 = 0.8e1 * t346 * t347 * t460 + 0.8e1 * t346 * t347 * t463 - 0.2e1 * t342 * t508;
  t576 = t175 * t32;
  t583 = t345 * t175;
  t589 = -0.2e1 * (t43 - t93 + t22 + t97 + t58 - t74 + 0.24e2) * t67 * t576 - 0.2e1 * (t43 + t93 + t22 + t97 - t58 - t74 + 0.24e2) * t67 * t176 + 0.4e1 * t583 * t326 - 0.4e1 * t583 * t325 * t34;
  t593 = (t179 - t58 - t24 + 0.8e1) * t206 * t32;
  t598 = (t179 + t58 - t24 + 0.8e1) * t206 * t10;
  t601 = t67 * t67;
  t603 = beta_y * t601 * t15;
  t609 = -0.4e1 * t34 * t347 * t603 - 0.2e1 * t342 * t593 - 0.2e1 * t342 * t598 + 0.4e1 * t348 * t603;
  t611 = (0.32e2 * t17 * t26 * t460 + 0.32e2 * t17 * t26 * t463 - 0.32e2 * t392 * t7 * t9) * t56 + (0.32e2 * t392 * t62 * t64 - 0.32e2 * t460 * t76 * t77 - 0.32e2 * t463 * t76 * t77) * t111 + (-0.32e2 * t114 * t116 * t392 - 0.32e2 * t120 * t121 * t460 - 0.32e2 * t120 * t121 * t463) * t140 + (0.32e2 * t143 * t144 * t392 - 0.64e2 * t144 * t460 * t58 - 0.64e2 * t144 * t463 * t58) * t168 + (-0.4e1 * t174 * t414 - 0.8e1 * t183 * t500 - 0.8e1 * t183 * t503) * t200 + (0.16e2 * t211 * t212 * t460 + 0.16e2 * t211 * t212 * t463 - 0.4e1 * t203 * t508) * t233 + (0.4e1 * t236 * t414 - 0.8e1 * t240 * t500 - 0.8e1 * t240 * t503) * t253 + (-0.16e2 * t259 * t260 * t460 - 0.16e2 * t259 * t260 * t463 + 0.4e1 * t256 * t508) * t276 + (-0.4e1 * t279 * t537 + 0.16e2 * t285 * t460 + 0.16e2 * t285 * t463) * t298 + (0.4e1 * t300 * t537 - 0.16e2 * t304 * t460 - 0.16e2 * t304 * t463) * t317 + t562 * t338 + t572 * t365 + t589 * t375 + t609 * t383;
  t614 = t19 * t10;
  t615 = t614 * t32;
  t618 = t389 * t16;
  t625 = t16 * t10;
  t629 = t397 * t16;
  t637 = t403 * t19;
  t643 = t10 * t32;
  t646 = t409 * t67;
  t652 = t176 * t32;
  t663 = t418 * t67;
  t693 = t435 * t340;
  t697 = -0.4e1 * t436 * t643 + 0.2e1 * t49 * t693 + 0.2e1 * t693;
  t701 = t601 * t175;
  t707 = 0.6e1 * t356 * t49 * t701 - 0.12e2 * t342 * t652 + 0.6e1 * t356 * t701;
  t709 = (-0.32e2 * t19 * t618 + 0.64e2 * t391 * t615 - 0.32e2 * t50 * t618) * t56 + (-0.64e2 * t32 * t398 * t625 + 0.32e2 * t49 * t629 + 0.32e2 * t629) * t111 + (-0.64e2 * t404 * t615 + 0.32e2 * t49 * t637 + 0.32e2 * t637) * t140 + (-0.8e1 * t410 * t643 + 0.4e1 * t49 * t646 + 0.4e1 * t646) * t200 + (-0.24e2 * t203 * t652 + 0.12e2 * t221 * t358 + 0.12e2 * t221 * t359) * t233 + (0.8e1 * t419 * t643 - 0.4e1 * t49 * t663 - 0.4e1 * t663) * t253 + (0.24e2 * t256 * t652 - 0.12e2 * t269 * t358 - 0.12e2 * t269 * t359) * t276 + (0.12e2 * t255 * t49 - 0.24e2 * t426 * t643 + 0.12e2 * t255) * t298 + (-0.12e2 * t202 * t49 + 0.24e2 * t430 * t643 - 0.12e2 * t202) * t317 + t697 * t338 + t707 * t365;
  t740 = t175 * t367;
  t742 = t175 * t370;
  t782 = (0.32e2 * t17 * t26 * t367 + 0.32e2 * t17 * t26 * t370) * t56 + (-0.32e2 * t367 * t76 * t77 - 0.32e2 * t370 * t76 * t77) * t111 + (-0.32e2 * t120 * t121 * t367 - 0.32e2 * t120 * t121 * t370) * t140 + (-0.64e2 * t144 * t367 * t58 - 0.64e2 * t144 * t370 * t58) * t168 + (-0.8e1 * t183 * t740 - 0.8e1 * t183 * t742) * t200 + (0.16e2 * t211 * t212 * t367 + 0.16e2 * t211 * t212 * t370) * t233 + (-0.8e1 * t240 * t740 - 0.8e1 * t240 * t742) * t253 + (-0.16e2 * t259 * t260 * t367 - 0.16e2 * t259 * t260 * t370) * t276 + (0.16e2 * t285 * t367 + 0.16e2 * t285 * t370) * t298 + (-0.16e2 * t304 * t367 - 0.16e2 * t304 * t370) * t317 - t374 * t338 - t382 * t365 + t337 * t375 + t364 * t383;
  t785 = t16 * t19;
  t810 = t224 * t32;
  t812 = t224 * t10;
  t847 = (-0.64e2 * t10 * t390 * t785 - 0.64e2 * t32 * t390 * t785) * t56 + (0.64e2 * t16 * t32 * t398 + 0.64e2 * t398 * t625) * t111 + (0.64e2 * t19 * t32 * t404 + 0.64e2 * t404 * t614) * t140 + (0.8e1 * t10 * t410 + 0.8e1 * t32 * t410) * t200 + (0.24e2 * t221 * t810 + 0.24e2 * t221 * t812) * t233 + (-0.8e1 * t10 * t419 - 0.8e1 * t32 * t419) * t253 + (-0.24e2 * t269 * t810 - 0.24e2 * t269 * t812) * t276 + (0.24e2 * t10 * t426 + 0.24e2 * t32 * t426) * t298 + (-0.24e2 * t10 * t430 - 0.24e2 * t32 * t430) * t317 - t446 * t338 - t453 * t365 + 0.4e1 * t436 * t392 * t375 + 0.12e2 * t342 * t414 * t383;
  t851 = 0.8e1 * t88 * t90;
  t861 = beta_y * t16;
  t869 = 0.13e2 * t21;
  t870 = 0.11e2 * t58;
  t871 = 0.4e1 * beta_y;
  t872 = -t91 + t94 + t69 - t869 + t97 + t70 + t870 - t67 - t871 - t24 + 0.4e1;
  t877 = 0.18e2 * t58;
  t878 = 0.13e2 * s;
  t883 = t75 * t16;
  t896 = -t91 + t94 - t69 + t869 - t97 + t70 + t870 + t67 - t871 + t24 - 0.4e1;
  t900 = t119 * t19;
  t925 = t39 * t67 * t41;
  t926 = 0.4e1 * t93;
  t927 = 0.3e1 * t69;
  t928 = 0.8e1 * t58;
  t929 = 0.6e1 * s;
  t934 = 0.3e1 * t91;
  t935 = 0.6e1 * t21;
  t936 = 0.3e1 * t70;
  t937 = 0.20e2 * t58;
  t938 = 0.10e2 * s;
  t943 = t182 * t67;
  t969 = t239 * t67;
  t980 = t345 * t18;
  t988 = 0.2e1 * t93;
  t997 = t284 * s;
  t1013 = t303 * s;
  t1026 = (-0.8e1 * (-t170 + t851 + t388 + t171 - t101 + 0.12e2) * t16 * t9 * t32 - 0.8e1 * (-t170 + t851 - t388 + t171 - t101 - 0.12e2) * t6 * t614 - 0.16e2 * t861 * t28 + 0.16e2 * t861 * t26 * t34) * t56 + (0.8e1 * t872 * t16 * t64 * t32 + 0.8e1 * (-t69 + t46 + t71 - t877 - t67 + t172 + t878 - 0.20e2) * t6 * t65 + 0.16e2 * t883 * t64 * t27 - 0.16e2 * t883 * t84) * t111 + (0.8e1 * (-t69 + t46 - t71 + t877 - t67 - t172 + t878 - 0.20e2) * t9 * t116 * t32 + 0.8e1 * t896 * t19 * t117 + 0.16e2 * t900 * t116 * t27 - 0.16e2 * t900 * t128) * t140 + (-0.16e2 * (t21 - t58 - t60 + 0.4e1) * s * t144 * t32 - 0.16e2 * (t21 + t58 - t60 + 0.4e1) * s * t145 + 0.32e2 * t70 * t148 - 0.32e2 * t70 * t144 * t34) * t168 + (0.2e1 * (t925 - t91 + t926 - t927 + t179 + t70 - t928 + t100 + t101 - t929) * s * t576 + 0.2e1 * (t925 - t934 + t926 + t69 + t935 + t936 - t937 - t100 + t101 + t938) * s * t176 + 0.4e1 * t943 * t184 - 0.4e1 * t943 * t189) * t200 + (0.8e1 * t212 * t34 * t346 + 0.4e1 * t203 * t593 + 0.4e1 * t203 * t598 - 0.8e1 * t213 * t346) * t233 + (0.2e1 * (t925 + t934 - t926 + t69 + t935 - t936 + t937 - t100 - t101 + t938) * s * t576 + 0.2e1 * (t925 + t91 - t926 - t927 + t179 - t70 + t928 + t100 - t101 - t929) * s * t176 + 0.4e1 * t969 * t184 - 0.4e1 * t969 * t189) * t253 + (-0.8e1 * t260 * t34 * t980 - 0.4e1 * t256 * t593 - 0.4e1 * t256 * t598 + 0.8e1 * t261 * t980) * t276 + (0.4e1 * (t988 - t46 - t59 + t172 + t74 - 0.24e2) * s * t576 + 0.4e1 * t5 * (t179 - t58 - s + 0.8e1) * t281 - 0.8e1 * t997 * t184 + 0.8e1 * t997 * t189) * t298 + (-0.4e1 * t8 * (t179 + t58 - s + 0.8e1) * t280 * t32 - 0.4e1 * (t988 + t46 - t59 + t172 - t74 + 0.24e2) * s * t176 + 0.8e1 * t1013 * t184 - 0.8e1 * t1013 * t189) * t317 - t589 * t338 - t609 * t365 + t562 * t375 + t572 * t383;
  return(At_fH_re * t385 * PREF_V_CA + At_fA_re * t455 * PREF_V_CA + Bt_fH_re * t611 * PREF_V_CA + Bt_fA_re * t709 * PREF_V_CA + At_fH_im * t782 * PREF_V_CA + At_fA_im * t847 * PREF_V_CA + Bt_fH_im * t1026 * PREF_V_CA + Bt_fA_im * (t375 * t697 + t383 * t707) * PREF_V_CA);
}
