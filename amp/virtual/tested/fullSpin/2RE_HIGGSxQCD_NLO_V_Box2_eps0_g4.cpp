double Eval_V_B2 (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  c_double& cg  = I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_0;
  c_double& cg1 = I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_1;
  c_double& cg3 = I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_0;
  c_double& cg5 = I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_1;

  double t1;
  double t10;
  double t1000;
  double t1001;
  double t1016;
  double t1031;
  double t1037;
  double t1038;
  double t1077;
  double t1094;
  double t11;
  double t112;
  double t1128;
  double t114;
  double t1146;
  double t12;
  double t122;
  double t124;
  double t128;
  double t129;
  double t13;
  double t134;
  double t136;
  double t139;
  double t145;
  double t146;
  double t148;
  double t149;
  double t150;
  double t151;
  double t154;
  double t155;
  double t157;
  double t158;
  double t159;
  double t16;
  double t161;
  double t162;
  double t163;
  double t164;
  double t165;
  double t166;
  double t167;
  double t168;
  double t17;
  double t170;
  double t171;
  double t172;
  double t178;
  double t182;
  double t183;
  double t184;
  double t185;
  double t186;
  double t187;
  double t188;
  double t189;
  double t19;
  double t190;
  double t191;
  double t192;
  double t199;
  double t2;
  double t20;
  double t201;
  double t202;
  double t203;
  double t208;
  double t209;
  double t213;
  double t214;
  double t215;
  double t224;
  double t225;
  double t226;
  double t233;
  double t235;
  double t236;
  double t24;
  double t240;
  double t241;
  double t242;
  double t243;
  double t244;
  double t248;
  double t255;
  double t26;
  double t260;
  double t266;
  double t267;
  double t269;
  double t27;
  double t271;
  double t272;
  double t273;
  double t274;
  double t279;
  double t283;
  double t290;
  double t292;
  double t293;
  double t294;
  double t297;
  double t299;
  double t301;
  double t303;
  double t304;
  double t309;
  double t312;
  double t313;
  double t314;
  double t317;
  double t321;
  double t323;
  double t324;
  double t329;
  double t33;
  double t333;
  double t34;
  double t340;
  double t342;
  double t345;
  double t347;
  double t349;
  double t351;
  double t358;
  double t359;
  double t36;
  double t365;
  double t368;
  double t37;
  double t376;
  double t378;
  double t384;
  double t385;
  double t387;
  double t390;
  double t398;
  double t399;
  double t401;
  double t402;
  double t403;
  double t407;
  double t408;
  double t410;
  double t411;
  double t415;
  double t419;
  double t422;
  double t423;
  double t426;
  double t427;
  double t429;
  double t431;
  double t436;
  double t446;
  double t447;
  double t449;
  double t45;
  double t452;
  double t456;
  double t457;
  double t458;
  double t460;
  double t461;
  double t462;
  double t465;
  double t47;
  double t471;
  double t474;
  double t475;
  double t477;
  double t478;
  double t48;
  double t483;
  double t484;
  double t486;
  double t489;
  double t495;
  double t499;
  double t5;
  double t500;
  double t502;
  double t51;
  double t511;
  double t515;
  double t516;
  double t518;
  c_double t520;
  double t521;
  double t523;
  double t524;
  double t53;
  double t531;
  double t54;
  double t542;
  double t545;
  double t55;
  double t56;
  double t563;
  double t57;
  double t577;
  double t58;
  double t580;
  double t59;
  double t598;
  double t6;
  double t60;
  double t600;
  double t608;
  double t61;
  double t611;
  double t62;
  double t63;
  double t64;
  double t641;
  double t659;
  double t664;
  double t668;
  double t69;
  double t70;
  double t701;
  double t705;
  double t707;
  double t72;
  double t73;
  double t74;
  double t747;
  double t750;
  double t753;
  double t756;
  double t764;
  double t767;
  double t775;
  double t777;
  double t778;
  double t79;
  double t790;
  double t8;
  double t80;
  double t81;
  double t814;
  double t82;
  double t826;
  double t828;
  double t83;
  double t836;
  double t838;
  double t841;
  double t844;
  double t85;
  double t850;
  double t851;
  double t856;
  double t86;
  double t87;
  double t872;
  double t884;
  double t887;
  double t889;
  double t894;
  double t90;
  double t903;
  double t904;
  double t91;
  double t910;
  double t924;
  double t929;
  double t932;
  double t934;
  double t935;
  double t938;
  double t939;
  double t942;
  double t946;
  double t949;
  double t952;
  double t957;
  double t961;
  double t965;
  double t970;
  double t972;
  double t997;
  double t998;
  double t999;
  t1 = beta_y * beta_y;
  t2 = t1 * s;
  t5 = (t2 - s + 0.2e1) * FH0 * s;
  t6 = beta_y * s;
  t8 = 0.1e1 / (t6 - s + 0.2e1);
  t10 = 0.1e1 / (t6 + s - 0.2e1);
  t11 = t8 * t10;
  t12 = sp(s1, k2);
  t13 = t11 * t12;
  t16 = sp(s1, p1);
  t17 = t16 * beta_y;
  t19 = s * t8;
  t20 = t19 * t10;
  t24 = sp(k1, s2);
  t26 = sp(s2, p1);
  t27 = t26 * beta_y;
  t33 = 0.2e1 - s;
  t34 = t2 - s + 0.4e1;
  t36 = t33 * t34 * FH0;
  t37 = sp(s1, s2);
  t45 = RE(I1_MT2_MU2_0);
  t47 = t1 * t1;
  t48 = t47 * s;
  t51 = 0.2e1 * s;
  t53 = FA0 * (t48 - 0.3e1 * t2 - 0.4e1 * t1 + t51 - 0.4e1);
  t54 = beta_y + 0.1e1;
  t55 = t54 * t54;
  t56 = 0.1e1 / t55;
  t57 = t53 * t56;
  t58 = beta_y - 0.1e1;
  t59 = t58 * t58;
  t60 = 0.1e1 / t59;
  t61 = 0.1e1 / s;
  t62 = t60 * t61;
  t63 = EPS_(k1, k2, s1, s2);
  t64 = t62 * t63;
  t69 = (0.2e1 * t1 - s + 0.2e1) * FH0;
  t70 = 0.1e1 / t58;
  t72 = 0.1e1 / t54;
  t73 = 0.4e1 - s;
  t74 = 0.1e1 / t73;
  t79 = s * s;
  t80 = t47 * t79;
  t81 = t1 * t79;
  t82 = 0.4e1 * t81;
  t83 = 0.24e2 * t2;
  t85 = 0.3e1 * t79;
  t86 = 0.8e1 * s;
  t87 = t80 - t82 + t83 - 0.64e2 * t1 + t85 - t86;
  t90 = beta_y * t87 * FH0 * t74;
  t91 = t56 * t60;
  t112 = (t1 * t47 * t79 + 0.24e2 * s - 0.32e2 * t47 + 0.16e2 * t48 - 0.4e1 * t79 - 0.4e1 * t80 + 0.7e1 * t81 - t83 - 0.32e2) * FH0;
  t114 = t91 * t37;
  t122 = RE(I2_MT2_0_MT2_MU2_0);
  t124 = s * FH0;
  t128 = t17 * t74;
  t129 = t124 * t128;
  t134 = beta_y * t74;
  t136 = t124 * t26 * t134 * t12;
  t139 = (t2 - t51 + 0.8e1) * t74;
  t145 = (-0.64e2 * t12 * t124 * t74 - 0.32e2 * t129) * t24 + 0.32e2 * t136 - 0.16e2 * t124 * t139 * t37 + 0.16e2 * t124 * t139;
  t146 = RE(I2_S12_0_0_MU2_0);
  t148 = 0.4e1 * t6;
  t149 = 0.3e1 * s;
  t150 = t2 - t148 + t149 - 0.4e1;
  t151 = t150 * FA0;
  t154 = 0.2e1 * t6;
  t155 = 0.4e1 * beta_y;
  t157 = (t2 + t154 + t155 - t149 + 0.4e1) * FH0;
  t158 = t8 * t70;
  t159 = t158 * t12;
  t161 = t1 * beta_y;
  t162 = t161 * t79;
  t163 = 0.3e1 * t81;
  t164 = 0.16e2 * t2;
  t165 = beta_y * t79;
  t166 = 0.3e1 * t165;
  t167 = 0.24e2 * t6;
  t168 = 0.32e2 * beta_y;
  t170 = FH0 * (t162 - t163 + t164 + t166 - t167 - t79 + t168 + t86 - 0.16e2);
  t171 = t170 * t60;
  t172 = t8 * t61;
  t178 = t26 * t12;
  t182 = 0.2e1 * t162;
  t183 = t161 * s;
  t184 = 0.12e2 * t183;
  t185 = 0.6e1 * t165;
  t186 = 0.16e2 * t1;
  t187 = 0.28e2 * t6;
  t188 = 0.16e2 * beta_y;
  t189 = 0.16e2 * s;
  t190 = t80 - t182 + t184 + t82 - t83 - t185 + t186 + t187 + t85 - t188 - t189 + 0.16e2;
  t191 = FH0 * t190;
  t192 = t60 * t8;
  t199 = RE(I2_T11_0_MT2_MU2_0);
  t201 = t2 + t148 + t149 - 0.4e1;
  t202 = t201 * FA0;
  t203 = t61 * t56;
  t208 = (t2 - t154 - t155 - t149 + 0.4e1) * FH0;
  t209 = t10 * t72;
  t213 = FH0 * (t162 + t163 - t164 + t166 - t167 + t79 + t168 - t86 + 0.16e2);
  t214 = t213 * t56;
  t215 = t10 * t61;
  t224 = t80 + t182 - t184 + t82 - t83 + t185 + t186 - t187 + t85 + t188 - t189 + 0.16e2;
  t225 = FH0 * t224;
  t226 = t56 * t10;
  t233 = RE(I2_T12_0_MT2_MU2_0);
  t235 = FA0 * t79;
  t236 = 0.1e1 / t34;
  t240 = t33 * FH0;
  t241 = t34 * t34;
  t242 = 0.1e1 / t241;
  t243 = t79 * t242;
  t244 = t243 * t12;
  t248 = FH0 * t79;
  t255 = t33 * t26 * beta_y;
  t260 = t33 * t236;
  t266 = 0.256e3 * t235 * t236 * t63 + (-0.128e3 * beta_y * t16 * t242 * t248 * t33 - 0.128e3 * t240 * t244) * t24 + 0.128e3 * t255 * t248 * t242 * t12 - 0.64e2 * t260 * t248 * t37 + 0.64e2 * t260 * t248;
  t267 = RE(I3_0_0_S12_0_0_0_MU2_1);
  t269 = s * FA0;
  t271 = 0.64e2 * t269 * t63;
  t272 = t124 * t12;
  t273 = beta_y - 0.2e1;
  t274 = t273 * t16;
  t279 = t273 * t26;
  t283 = t2 - t154 - s + 0.4e1;
  t290 = RE(I3_0_T11_MT2_0_0_MT2_MU2_0);
  t292 = t59 * FA0;
  t293 = t79 * t236;
  t294 = t293 * t63;
  t297 = t2 + t6 - t51 + 0.4e1;
  t299 = t59 * t297 * FH0;
  t301 = t6 - s + 0.4e1;
  t303 = t59 * t301 * FH0;
  t304 = t243 * t16;
  t309 = t243 * t178;
  t312 = t59 * t33;
  t313 = t312 * FH0;
  t314 = t293 * t37;
  t317 = t248 * t236;
  t321 = RE(I3_0_T11_MT2_0_0_MT2_MU2_1);
  t323 = beta_y + 0.2e1;
  t324 = t323 * t16;
  t329 = t323 * t26;
  t333 = t2 + t154 - s + 0.4e1;
  t340 = RE(I3_0_T12_MT2_0_0_MT2_MU2_0);
  t342 = t55 * FA0;
  t345 = t2 - t6 - t51 + 0.4e1;
  t347 = t55 * t345 * FH0;
  t349 = t6 + s - 0.4e1;
  t351 = t55 * t349 * FH0;
  t358 = t55 * t33;
  t359 = t358 * FH0;
  t365 = RE(I3_0_T12_MT2_0_0_MT2_MU2_1);
  t368 = (0.4e1 - t149) * t74;
  t376 = 0.4e1 * t2;
  t378 = (t376 - t85 + t189 - 0.16e2) * t74;
  t384 = (0.16e2 * t12 * t124 * t368 - 0.64e2 * t129) * t24 + 0.64e2 * t136 - 0.8e1 * t124 * t378 * t37 + 0.8e1 * t124 * t378;
  t385 = RE(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t387 = ((0.64e2 * FH0 * t17 * t20 - 0.32e2 * t13 * t5) * t24 - 0.64e2 * t27 * FH0 * t19 * t10 * t12 - 0.16e2 * t36 * t19 * t10 * t37 + 0.16e2 * t36 * t20) * t45 + (0.128e3 * t57 * t64 + (0.128e3 * t12 * t69 * t70 * t72 * t74 + 0.32e2 * t16 * t61 * t90 * t91) * t24 - 0.32e2 * t90 * t91 * t61 * t26 * t12 + 0.16e2 * t112 * t74 * t114 - 0.16e2 * t112 * t74 * t56 * t60) * t122 + t145 * t146 + (-0.64e2 * t151 * t64 + (0.16e2 * t16 * t171 * t172 - 0.16e2 * t157 * t159) * t24 - 0.16e2 * t171 * t172 * t178 + 0.8e1 * t191 * t192 * t37 - 0.8e1 * t191 * t192) * t199 + (-0.64e2 * t202 * t203 * t63 + (-0.16e2 * t12 * t208 * t209 - 0.16e2 * t16 * t214 * t215) * t24 + 0.16e2 * t214 * t215 * t178 - 0.8e1 * t225 * t226 * t37 + 0.8e1 * t225 * t226) * t233 + t266 * t267 + (-t271 + (0.16e2 * t124 * t274 + 0.16e2 * t272) * t24 - 0.16e2 * t124 * t279 * t12 + 0.8e1 * t124 * t283 * t37 - 0.8e1 * t124 * t283) * t290 + (-0.128e3 * t292 * t294 + (0.32e2 * t244 * t299 - 0.32e2 * t303 * t304) * t24 + 0.32e2 * t303 * t309 + 0.32e2 * t313 * t314 - 0.32e2 * t312 * t317) * t321 + (-t271 + (0.16e2 * t124 * t324 + 0.16e2 * t272) * t24 - 0.16e2 * t124 * t329 * t12 + 0.8e1 * t124 * t333 * t37 - 0.8e1 * t124 * t333) * t340 + (-0.128e3 * t342 * t294 + (0.32e2 * t244 * t347 - 0.32e2 * t304 * t351) * t24 + 0.32e2 * t351 * t309 + 0.32e2 * t359 * t314 - 0.32e2 * t358 * t317) * t365 + t384 * t385;
  t390 = t12 * t24;
  t398 = -0.8e1 * t248 * t37 * t73 + 0.32e2 * t235 * t63 - 0.16e2 * t248 * t390 + 0.8e1 * t248 * t73;
  t399 = RE(cg);
  t401 = t79 * s;
  t402 = FA0 * t401;
  t403 = t59 * t236;
  t407 = t401 * t242;
  t408 = t407 * t12;
  t410 = t407 * t16;
  t411 = t303 * t410;
  t415 = t407 * t178;
  t419 = t401 * t236 * t37;
  t422 = FH0 * t401;
  t423 = t422 * t236;
  t426 = 0.64e2 * t402 * t403 * t63 + (-0.16e2 * t299 * t408 + 0.16e2 * t411) * t24 - 0.16e2 * t303 * t415 - 0.16e2 * t313 * t419 + 0.16e2 * t312 * t423;
  t427 = RE(cg1);
  t429 = RE(cg3);
  t431 = t55 * t236;
  t436 = t351 * t410;
  t446 = 0.64e2 * t402 * t431 * t63 + (-0.16e2 * t347 * t408 + 0.16e2 * t436) * t24 - 0.16e2 * t351 * t415 - 0.16e2 * t359 * t419 + 0.16e2 * t358 * t423;
  t447 = RE(cg5);
  t449 = EPS_(k1, k2, p1, s2);
  t452 = EPS_(k1, k2, p1, s1);
  t456 = beta_y * t124 * t452 * t74 + t124 * t134 * t449;
  t457 = 0.32e2 * t456;
  t458 = IM(I2_S12_0_0_MU2_0);
  t460 = t33 * beta_y;
  t461 = t460 * FH0;
  t462 = t243 * t452;
  t465 = t243 * t449;
  t471 = t12 * t236;
  t474 = -0.256e3 * t235 * t236 * t24 - 0.256e3 * t235 * t471 + 0.128e3 * t461 * t462 + 0.128e3 * t461 * t465;
  t475 = IM(I3_0_0_S12_0_0_0_MU2_1);
  t477 = 0.64e2 * t456;
  t478 = IM(I3_S12_MT2_MT2_0_0_MT2_MU2_0);
  t483 = -0.32e2 * t12 * t235 - 0.32e2 * t235 * t24;
  t484 = IM(cg);
  t486 = t407 * t452;
  t489 = t407 * t449;
  t495 = t59 * t12;
  t499 = -0.64e2 * t236 * t402 * t495 - 0.64e2 * t24 * t402 * t403 - 0.16e2 * t303 * t486 - 0.16e2 * t303 * t489;
  t500 = IM(cg1);
  t502 = IM(cg3);
  t511 = t55 * t12;
  t515 = -0.64e2 * t236 * t402 * t511 - 0.64e2 * t24 * t402 * t431 - 0.16e2 * t351 * t486 - 0.16e2 * t351 * t489;
  t516 = IM(cg5);
  t518 = t398 * t399 + t398 * t429 + t426 * t427 + t446 * t447 + t457 * t458 + t474 * t475 + t477 * t478 + t483 * t484 + t483 * t502 + t499 * t500 + t515 * t516;
  t520 = DenS(s, mH, GammaH);
  t521 = RE(t520);
  t523 = beta_y * FH0;
  t524 = t523 * s;
  t531 = IM(t520);
  t542 = t62 * t24;
  t545 = t53 * t12;
  t563 = t150 * t12;
  t577 = t203 * t24;
  t580 = t201 * t12;
  t598 = 0.64e2 * t269 * t24;
  t600 = 0.64e2 * t269 * t12;
  t608 = t293 * t24;
  t611 = t235 * t236;
  t641 = -t483 * t531;
  t659 = (t387 + t518) * t521 + (0.64e2 * t11 * t449 * t524 + 0.64e2 * t11 * t452 * t524) * t531 * t45 + (0.32e2 * t449 * t61 * t90 * t91 + 0.32e2 * t452 * t61 * t90 * t91 + 0.128e3 * t545 * t61 * t91 + 0.128e3 * t542 * t57) * t531 * t122 - t457 * t531 * t146 + (-0.64e2 * FA0 * t563 * t60 * t61 + 0.16e2 * t171 * t172 * t449 + 0.16e2 * t171 * t172 * t452 - 0.64e2 * t151 * t542) * t531 * t199 + (-0.64e2 * FA0 * t56 * t580 * t61 - 0.16e2 * t214 * t215 * t449 - 0.16e2 * t214 * t215 * t452 - 0.64e2 * t202 * t577) * t531 * t233 - t474 * t531 * t267 + (0.16e2 * t124 * t273 * t449 + 0.16e2 * t124 * t273 * t452 - t598 - t600) * t531 * t290 + (-0.128e3 * t292 * t608 - 0.32e2 * t303 * t462 - 0.32e2 * t303 * t465 - 0.128e3 * t495 * t611) * t531 * t321 + (0.16e2 * t124 * t323 * t449 + 0.16e2 * t124 * t323 * t452 - t598 - t600) * t531 * t340 + (-0.128e3 * t342 * t608 - 0.32e2 * t351 * t462 - 0.32e2 * t351 * t465 - 0.128e3 * t511 * t611) * t531 * t365 - t477 * t531 * t385 + t641 * t399 - t499 * t531 * t427 + t641 * t429 - t515 * t531 * t447 + (t145 * t458 + t266 * t475 + t384 * t478 + t398 * t484 + t398 * t502 + t426 * t500 + t446 * t516) * t531;
  t664 = EPS_(k1, p1, s1, s2);
  t668 = EPS_(k2, p1, s1, s2);
  t701 = t124 * t134 * t664;
  t705 = t124 * t668 * beta_y * t74;
  t707 = -0.64e2 * t124 * t63 * t74 - 0.32e2 * t701 - 0.32e2 * t705;
  t747 = t243 * t63;
  t750 = t243 * t664;
  t753 = t243 * t668;
  t756 = t471 * t24;
  t764 = 0.128e3 * t236 * t37 * t402 - 0.256e3 * t235 * t756 + 0.128e3 * t236 * t402 - 0.128e3 * t240 * t747 - 0.128e3 * t461 * t750 - 0.128e3 * t461 * t753;
  t767 = 0.16e2 * t124 * t63;
  t775 = 0.64e2 * t269 * t390;
  t777 = 0.32e2 * t235 * t37;
  t778 = 0.32e2 * t235;
  t790 = t403 * t37;
  t814 = t431 * t37;
  t826 = 0.16e2 * t124 * t368 * t63 - 0.64e2 * t701 - 0.64e2 * t705;
  t828 = (-0.32e2 * t11 * t5 * t63 + 0.64e2 * t11 * t524 * t664 + 0.64e2 * t11 * t524 * t668) * t45 + (0.128e3 * t63 * t69 * t70 * t72 * t74 - 0.128e3 * t24 * t545 * t61 * t91 + 0.32e2 * t61 * t664 * t90 * t91 + 0.32e2 * t61 * t668 * t90 * t91 + 0.64e2 * t114 * t53 + 0.64e2 * t53 * t91) * t122 + t707 * t146 + (0.64e2 * FA0 * t542 * t563 - 0.32e2 * t151 * t37 * t60 - 0.16e2 * t157 * t158 * t63 + 0.16e2 * t171 * t172 * t664 + 0.16e2 * t171 * t172 * t668 - 0.32e2 * t151 * t60) * t199 + (0.64e2 * FA0 * t577 * t580 - 0.32e2 * t202 * t37 * t56 - 0.16e2 * t208 * t209 * t63 - 0.16e2 * t214 * t215 * t664 - 0.16e2 * t214 * t215 * t668 - 0.32e2 * t202 * t56) * t233 + t764 * t267 + (0.16e2 * t124 * t273 * t664 + 0.16e2 * t124 * t273 * t668 + t767 + t775 - t777 - t778) * t290 + (0.128e3 * FA0 * t495 * t608 + 0.32e2 * t299 * t747 - 0.32e2 * t303 * t750 - 0.32e2 * t303 * t753 - 0.64e2 * t402 * t403 - 0.64e2 * t402 * t790) * t321 + (0.16e2 * t124 * t323 * t664 + 0.16e2 * t124 * t323 * t668 + t767 + t775 - t777 - t778) * t340 + (0.128e3 * FA0 * t511 * t608 + 0.32e2 * t347 * t747 - 0.32e2 * t351 * t750 - 0.32e2 * t351 * t753 - 0.64e2 * t402 * t431 - 0.64e2 * t402 * t814) * t365 + t826 * t385;
  t836 = -0.32e2 * t235 * t390 - 0.16e2 * t248 * t63 + 0.16e2 * t37 * t402 + 0.16e2 * t402;
  t838 = t407 * t63;
  t841 = t407 * t664;
  t844 = t407 * t668;
  t850 = t79 * t79;
  t851 = t850 * FA0;
  t856 = -0.64e2 * t402 * t59 * t756 - 0.16e2 * t299 * t838 + 0.16e2 * t303 * t841 + 0.16e2 * t303 * t844 + 0.32e2 * t403 * t851 + 0.32e2 * t790 * t851;
  t872 = -0.64e2 * t402 * t55 * t756 - 0.16e2 * t347 * t838 + 0.16e2 * t351 * t841 + 0.16e2 * t351 * t844 + 0.32e2 * t431 * t851 + 0.32e2 * t814 * t851;
  t884 = t248 * t128;
  t887 = t248 * t27 * t74;
  t889 = -0.8e1 * t124 * (t2 + t6 + 0.8e1) * t74 * t24 - 0.8e1 * t124 * (t2 - t6 + 0.8e1) * t74 * t12 - 0.16e2 * t884 + 0.16e2 * t887;
  t894 = t243 * t24;
  t903 = FH0 * t242;
  t904 = t903 * t16;
  t910 = -0.32e2 * t33 * (t2 + t6 + 0.4e1) * FH0 * t894 - 0.32e2 * t33 * (t2 - t6 + 0.4e1) * FH0 * t244 - 0.64e2 * t460 * t401 * t904 + 0.64e2 * t255 * t422 * t242;
  t924 = -0.16e2 * t124 * (t2 + t6 + t149 - 0.4e1) * t74 * t24 - 0.16e2 * t124 * (t2 - t6 + t149 - 0.4e1) * t74 * t12 - 0.32e2 * t884 + 0.32e2 * t887;
  t929 = -0.16e2 * t12 * t248 - 0.16e2 * t24 * t248;
  t932 = -t81 + t376 + t79 - 0.12e2 * s + 0.16e2;
  t934 = t59 * t932 * FH0;
  t935 = t407 * t24;
  t938 = 0.2e1 * t165;
  t939 = 0.4e1 * s;
  t942 = t59 * (-t81 + t376 + t938 - t79 - t939 + 0.16e2) * FH0;
  t946 = t59 * t850 * t301;
  t949 = t903 * t26;
  t952 = -0.4e1 * t408 * t942 + 0.8e1 * t904 * t946 - 0.4e1 * t934 * t935 - 0.8e1 * t946 * t949;
  t957 = t55 * (-t81 + t376 - t938 - t79 - t939 + 0.16e2) * FH0;
  t961 = t55 * t932 * FH0;
  t965 = t55 * t850 * t349;
  t970 = -0.4e1 * t408 * t961 + 0.8e1 * t904 * t965 - 0.4e1 * t935 * t957 - 0.8e1 * t949 * t965;
  t972 = t399 * t836 + t427 * t856 + t429 * t836 + t447 * t872 + t458 * t889 + t475 * t910 + t478 * t924 + t484 * t929 + t500 * t952 + t502 * t929 + t516 * t970;
  t997 = t47 * beta_y * t79;
  t998 = 0.4e1 * t162;
  t999 = 0.24e2 * t183;
  t1000 = 0.32e2 * t161;
  t1001 = 0.32e2 * t1;
  t1016 = t523 * t87;
  t1031 = 0.12e2 * t2;
  t1037 = 0.32e2 * t6;
  t1038 = 0.20e2 * s;
  t1077 = 0.3e1 * t6;
  t1094 = t407 * t26;
  t1128 = -t929 * t531;
  t1146 = (t828 + t972) * t521 + (0.16e2 * FH0 * s * t11 * t24 * t345 + 0.32e2 * FH0 * t10 * t27 * t79 * t8 + 0.16e2 * FH0 * s * t13 * t297 - 0.32e2 * FH0 * t11 * t16 * t165) * t531 * t45 + (-0.8e1 * FH0 * (t997 - t998 + t999 - t1000 + t166 - t1001 - t167 + t168 + t189 - 0.32e2) * t74 * t72 * t60 * t24 - 0.8e1 * FH0 * (t997 - t998 + t999 - t1000 + t166 + t1001 - t167 + t168 - t189 + 0.32e2) * t74 * t56 * t70 * t12 - 0.16e2 * t1016 * t91 * t74 * t16 + 0.16e2 * t1016 * t91 * t74 * t26) * t531 * t122 - t889 * t531 * t146 + (-0.4e1 * FH0 * (t80 - t182 + t184 - t1031 + t938 + t186 + t148 - t79 + t188 - t939) * t192 * t24 - 0.4e1 * (t162 - t163 + t1031 + t166 - t1037 - t79 + t188 + t1038 - 0.32e2) * FH0 * t159 - 0.8e1 * t170 * t192 * t16 + 0.8e1 * t170 * t192 * t26) * t531 * t199 + (0.4e1 * (t162 + t163 - t1031 + t166 - t1037 + t79 + t188 - t1038 + 0.32e2) * FH0 * t209 * t24 + 0.4e1 * FH0 * (t80 + t182 - t184 - t1031 - t938 + t186 - t148 - t79 - t188 - t939) * t226 * t12 + 0.8e1 * t213 * t226 * t16 - 0.8e1 * t213 * t226 * t26) * t531 * t233 - t910 * t531 * t267 + (-0.4e1 * t124 * t345 * t24 - 0.4e1 * t124 * (t2 - t1077 + t51 + 0.4e1) * t12 - 0.8e1 * t248 * t274 + 0.8e1 * t248 * t279) * t531 * t290 + (-0.16e2 * t1094 * t303 - 0.8e1 * t244 * t942 - 0.8e1 * t894 * t934 + 0.16e2 * t411) * t531 * t321 + (-0.4e1 * t124 * (t2 + t1077 + t51 + 0.4e1) * t24 - 0.4e1 * t124 * t297 * t12 - 0.8e1 * t248 * t324 + 0.8e1 * t248 * t329) * t531 * t340 + (-0.16e2 * t1094 * t351 - 0.8e1 * t244 * t961 - 0.8e1 * t894 * t957 + 0.16e2 * t436) * t531 * t365 - t924 * t531 * t385 + t1128 * t399 - t952 * t531 * t427 + t1128 * t429 - t970 * t531 * t447 + (t458 * t707 + t475 * t764 + t478 * t826 + t484 * t836 + t500 * t856 + t502 * t836 + t516 * t872) * t531;
  return(t1146 * PREF_V_CA_Bt + t659 * PREF_V_CA_At);
}
