double Eval_V_SE (PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t10;
  double t102;
  double t107;
  double t11;
  double t111;
  double t115;
  double t12;
  double t121;
  double t126;
  double t13;
  double t132;
  double t14;
  double t141;
  double t143;
  double t145;
  double t147;
  double t149;
  double t15;
  double t150;
  double t155;
  double t16;
  double t160;
  double t163;
  double t17;
  double t170;
  double t171;
  double t172;
  double t175;
  double t177;
  double t178;
  double t18;
  double t180;
  double t182;
  double t183;
  double t189;
  double t19;
  double t191;
  double t193;
  double t194;
  double t195;
  double t197;
  double t198;
  double t199;
  double t2;
  double t20;
  double t200;
  double t21;
  double t214;
  double t22;
  double t221;
  double t226;
  double t23;
  double t230;
  double t232;
  double t233;
  double t234;
  double t235;
  double t236;
  double t237;
  double t238;
  double t24;
  double t240;
  double t241;
  double t242;
  double t243;
  double t247;
  double t248;
  double t25;
  double t251;
  double t252;
  double t254;
  double t255;
  double t256;
  double t26;
  double t264;
  double t265;
  double t266;
  double t267;
  double t268;
  double t269;
  double t27;
  double t270;
  double t271;
  double t272;
  double t273;
  double t274;
  double t275;
  double t276;
  double t277;
  double t278;
  double t279;
  double t28;
  double t280;
  double t281;
  double t286;
  double t287;
  double t29;
  double t291;
  double t294;
  double t295;
  double t296;
  double t3;
  double t300;
  double t301;
  double t305;
  double t306;
  double t307;
  double t308;
  double t316;
  double t317;
  double t319;
  double t32;
  double t322;
  double t323;
  double t327;
  double t329;
  double t330;
  double t331;
  double t333;
  double t334;
  double t337;
  double t34;
  double t345;
  double t35;
  double t353;
  double t357;
  double t362;
  double t365;
  double t378;
  double t38;
  double t381;
  double t386;
  double t39;
  double t393;
  double t396;
  double t40;
  double t402;
  double t409;
  double t421;
  double t422;
  double t426;
  double t44;
  double t45;
  double t453;
  double t46;
  double t47;
  double t48;
  double t49;
  double t493;
  double t495;
  double t496;
  double t497;
  double t498;
  double t499;
  double t5;
  double t500;
  double t501;
  double t502;
  double t503;
  double t504;
  double t505;
  double t506;
  double t507;
  double t514;
  double t532;
  double t533;
  double t534;
  double t535;
  double t536;
  double t541;
  double t542;
  double t543;
  double t544;
  double t545;
  double t546;
  double t55;
  double t559;
  double t564;
  double t59;
  double t6;
  double t60;
  double t61;
  double t62;
  double t70;
  double t8;
  double t81;
  double t82;
  double t84;
  double t85;
  double t88;
  double t89;
  double t9;
  c_double t93;
  double t94;
  double t96;
  double t97;
  t1 = Z2(mt2);
  t2 = beta_y * beta_y;
  t3 = t2 * t2;
  t5 = t1 * t3 * s;
  t6 = Z2M(mt2);
  t8 = t6 * t3 * s;
  t9 = t1 * t2;
  t10 = 0.4e1 * t9;
  t11 = t6 * t2;
  t12 = 0.4e1 * t11;
  t13 = t11 * s;
  t14 = 0.4e1 * t1;
  t15 = 0.4e1 * t6;
  t16 = t1 * s;
  t17 = t5 - t8 + t10 - t12 + t13 + t14 - t15 - t16;
  t18 = FA0 * t17;
  t19 = 0.1e1 / s;
  t20 = t18 * t19;
  t21 = beta_y + 0.1e1;
  t22 = t21 * t21;
  t23 = 0.1e1 / t22;
  t24 = beta_y - 0.1e1;
  t25 = t24 * t24;
  t26 = 0.1e1 / t25;
  t27 = t23 * t26;
  t28 = EPS_(k1, k2, s1, s2);
  t29 = t27 * t28;
  t32 = t9 * s;
  t34 = 0.2e1 * t16;
  t35 = t6 * s;
  t38 = FH0 * (t10 - t12 + 0.2e1 * t32 - t13 + t14 - t15 - t34 + t35) * t19;
  t39 = sp(s1, k2);
  t40 = t27 * t39;
  t44 = FH0 * (t1 - t6);
  t45 = t44 * beta_y;
  t46 = 0.1e1 / t24;
  t47 = 0.1e1 / t21;
  t48 = t46 * t47;
  t49 = sp(s1, p1);
  t55 = sp(k1, s2);
  t59 = t47 * t19;
  t60 = sp(s2, p1);
  t61 = t60 * t39;
  t62 = t59 * t61;
  t70 = s * s;
  t81 = 0.2e1 * t1 * t70 + t11 * t70 - t6 * t70 - 0.2e1 * t70 * t9 + 0.16e2 * t1 - 0.16e2 * t11 + 0.4e1 * t13 - 0.12e2 * t16 + 0.8e1 * t35 + 0.4e1 * t5 - 0.16e2 * t6 - 0.4e1 * t8 + 0.16e2 * t9;
  t82 = FH0 * t81;
  t84 = sp(s1, s2);
  t85 = t27 * t84;
  t88 = t19 * t23;
  t89 = t88 * t26;
  t93 = DenS(s, mH, GammaH);
  t94 = RE(t93);
  t96 = t46 * t19;
  t97 = EPS_(k1, k2, p1, s1);
  t102 = EPS_(k1, k2, p1, s2);
  t107 = t27 * t55;
  t111 = FA0 * t39 * t17;
  t115 = IM(t93);
  t121 = EPS_(k1, p1, s1, s2);
  t126 = EPS_(k2, p1, s1, s2);
  t132 = t88 * t26 * t55;
  t141 = t2 * beta_y;
  t143 = t1 * t141 * s;
  t145 = t6 * t141 * s;
  t147 = t1 * beta_y * s;
  t149 = t6 * beta_y * s;
  t150 = t5 - t8 + t143 + t10 - t145 - t12 + t32 - t147 + t14 + t149 - t15 - t34 + t35;
  t155 = t5 - t8 - t143 + t10 + t145 - t12 + t32 + t147 + t14 - t149 - t15 - t34 + t35;
  t160 = t48 * t49;
  t163 = t48 * t60;
  t170 = t2 * s;
  t171 = t170 - s + 0.4e1;
  t172 = FA0 * t171;
  t175 = -0.3e1 * t170 + 0.2e1 * t2 - s + 0.2e1;
  t177 = t172 * t175 * t19;
  t178 = beta_y * s;
  t180 = 0.1e1 / (t178 - s + 0.2e1);
  t182 = 0.1e1 / (t178 + s - 0.2e1);
  t183 = t180 * t182;
  t189 = FH0 * (t2 + 0.1e1);
  t191 = t189 * t19 * t40;
  t193 = FH0 * beta_y;
  t194 = t2 * t70;
  t195 = 0.8e1 * s;
  t197 = (t194 - t70 + t195 - 0.8e1) * t182;
  t198 = t193 * t197;
  t199 = t180 * t46;
  t200 = t59 * t49;
  t214 = t3 * t70;
  t221 = FH0 * t171 * (t214 - 0.3e1 * t194 + 0.12e2 * t170 - 0.8e1 * t2 - 0.2e1 * t70 + t195 - 0.8e1) * t19;
  t226 = t27 * t183;
  t230 = RE(I1_MT2_MU2_0);
  t232 = t141 * t70;
  t233 = 0.3e1 * t232;
  t234 = 0.7e1 * t194;
  t235 = 0.4e1 * t170;
  t236 = beta_y * t70;
  t237 = 0.5e1 * t236;
  t238 = 0.4e1 * t178;
  t240 = FA0 * (t233 - t234 + t235 + t237 + t238 - t70 - t195 + 0.16e2);
  t241 = t240 * t180;
  t242 = t26 * t19;
  t243 = t242 * t28;
  t247 = (t178 - s + 0.4e1) * FH0;
  t248 = t242 * t39;
  t251 = 0.2e1 * t236;
  t252 = 0.8e1 * t178;
  t254 = (t194 - t251 - t252 + t70 + t195 - 0.16e2) * FH0;
  t255 = t254 * t19;
  t256 = t199 * t49;
  t264 = t70 * s;
  t265 = t3 * t264;
  t266 = t141 * t264;
  t267 = 0.3e1 * t266;
  t268 = 0.8e1 * t232;
  t269 = t2 * t264;
  t270 = 0.5e1 * t269;
  t271 = 0.8e1 * t194;
  t272 = beta_y * t264;
  t273 = 0.5e1 * t272;
  t274 = 0.16e2 * t170;
  t275 = 0.20e2 * t236;
  t276 = 0.2e1 * t264;
  t277 = 0.32e2 * t178;
  t278 = 0.20e2 * t70;
  t279 = 0.64e2 * s;
  t280 = t265 - t267 - t268 + t270 + t271 - t273 - t274 + t275 + t276 - t277 - t278 + t279 - 0.64e2;
  t281 = FH0 * t280;
  t286 = t180 * t26;
  t287 = t286 * t19;
  t291 = RE(I2_T11_0_MT2_MU2_0);
  t294 = FA0 * (t233 + t234 - t235 + t237 + t238 + t70 + t195 - 0.16e2);
  t295 = t294 * t182;
  t296 = t88 * t28;
  t300 = (t178 + s - 0.4e1) * FH0;
  t301 = t88 * t39;
  t305 = (t194 + t251 + t252 + t70 + t195 - 0.16e2) * FH0;
  t306 = t305 * t19;
  t307 = t47 * t182;
  t308 = t307 * t49;
  t316 = t265 + t267 + t268 + t270 + t271 + t273 - t274 - t275 + t276 + t277 - t278 + t279 - 0.64e2;
  t317 = FH0 * t316;
  t319 = t88 * t84;
  t322 = t182 * t23;
  t323 = t322 * t19;
  t327 = RE(I2_T12_0_MT2_MU2_0);
  t329 = t3 * s;
  t330 = 0.2e1 * t170;
  t331 = 0.4e1 * t2;
  t333 = FA0 * (t329 - t330 + t331 + s + 0.4e1);
  t334 = t333 * t26;
  t337 = t193 * t46;
  t345 = FH0 * (t329 - t330 + t331 - s + 0.4e1);
  t353 = t59 * t97;
  t357 = t59 * t102;
  t362 = t27 * t183 * t55;
  t365 = t175 * t39;
  t378 = t88 * t55;
  t381 = t333 * t39;
  t386 = t254 * t180;
  t393 = t242 * t55;
  t396 = t240 * t39;
  t402 = t305 * t182;
  t409 = t294 * t39;
  t421 = 0.128e3 * t189 * t26 * t296;
  t422 = t59 * t121;
  t426 = t59 * t126;
  t453 = t19 * t55;
  t493 = t3 * t2 * t264;
  t495 = t3 * beta_y * t264;
  t496 = 0.2e1 * t265;
  t497 = 0.16e2 * t214;
  t498 = 0.2e1 * t266;
  t499 = 0.8e1 * t329;
  t500 = t141 * s;
  t501 = 0.8e1 * t500;
  t502 = 0.40e2 * t170;
  t503 = 0.8e1 * t236;
  t504 = 0.32e2 * t2;
  t505 = 0.8e1 * t70;
  t506 = 0.32e2 * s;
  t507 = t493 + t495 - t496 + t497 - t498 - t499 + t268 + t269 - t501 - t271 + t272 + t502 - t503 - t504 + t252 - t505 + t506 - 0.32e2;
  t514 = t493 - t495 - t496 + t497 + t498 - t499 - t268 + t269 + t501 - t271 - t272 + t502 + t503 - t504 - t252 - t505 + t506 - 0.32e2;
  t532 = 0.2e1 * t272;
  t533 = 0.24e2 * t236;
  t534 = 0.48e2 * t178;
  t535 = 0.16e2 * t70;
  t536 = t265 - t498 - t268 + t532 - t274 + t533 - t264 - t534 - t535 + t279 - 0.64e2;
  t541 = 0.4e1 * t266;
  t542 = 0.6e1 * t269;
  t543 = 0.16e2 * t194;
  t544 = 0.4e1 * t272;
  t545 = 0.16e2 * t178;
  t546 = t265 - t541 - t268 + t542 + t543 - t544 - t274 - t503 + t264 - t545 + t506 - 0.64e2;
  t559 = t265 + t541 + t268 + t542 + t543 + t544 - t274 + t503 + t264 + t545 + t506 - 0.64e2;
  t564 = t265 + t498 + t268 - t532 - t274 - t533 - t264 + t534 - t535 + t279 - 0.64e2;
  return(((0.256e3 * t20 * t29 + (-0.512e3 * t19 * t45 * t48 * t49 - 0.128e3 * t38 * t40) * t55 + 0.512e3 * t44 * beta_y * t46 * t62 - 0.64e2 * t82 * t19 * t85 + 0.64e2 * t82 * t89) * t94 + (-0.512e3 * t102 * t45 * t47 * t96 - 0.512e3 * t45 * t47 * t96 * t97 + 0.256e3 * t107 * t20 + 0.256e3 * t111 * t89) * t115) * PREF_B_At + ((-0.512e3 * t121 * t45 * t47 * t96 - 0.512e3 * t126 * t45 * t47 * t96 - 0.256e3 * t111 * t132 + 0.128e3 * t18 * t27 + 0.128e3 * t18 * t85 - 0.128e3 * t29 * t38) * t94 + (0.128e3 * FH0 * t107 * t150 * t19 + 0.128e3 * FH0 * t155 * t19 * t40 + 0.256e3 * t160 * t45 - 0.256e3 * t163 * t45) * t115) * PREF_B_Bt + (((0.128e3 * t177 * t27 * t183 * t28 + (0.64e2 * t198 * t199 * t200 + 0.128e3 * t191) * t55 - 0.64e2 * t193 * t197 * t180 * t48 * t19 * t60 * t39 + 0.32e2 * t221 * t27 * t183 * t84 - 0.32e2 * t221 * t226) * t230 + (-0.32e2 * t241 * t243 + (0.32e2 * t247 * t248 - 0.16e2 * t255 * t256) * t55 + 0.16e2 * t255 * t199 * t61 - 0.8e1 * t281 * t180 * t242 * t84 + 0.8e1 * t281 * t287) * t291 + (-0.32e2 * t295 * t296 + (-0.32e2 * t300 * t301 + 0.16e2 * t306 * t308) * t55 - 0.16e2 * t306 * t307 * t61 + 0.8e1 * t317 * t182 * t319 - 0.8e1 * t317 * t323) * t327 + 0.64e2 * t334 * t296 + (-0.128e3 * t200 * t337 - 0.128e3 * t191) * t55 + 0.128e3 * t337 * t62 - 0.64e2 * t345 * t26 * t319 + 0.64e2 * t345 * t89) * t94 + (0.128e3 * t172 * t182 * t286 * t365 * t88 + 0.64e2 * t198 * t199 * t353 + 0.64e2 * t198 * t199 * t357 + 0.128e3 * t177 * t362) * t115 * t230 + (0.64e2 * t334 * t378 - 0.128e3 * t337 * t353 - 0.128e3 * t337 * t357 + 0.64e2 * t381 * t89) * t115 + (-0.16e2 * t102 * t386 * t96 - 0.16e2 * t386 * t96 * t97 - 0.32e2 * t241 * t393 - 0.32e2 * t287 * t396) * t115 * t291 + (-0.32e2 * t295 * t378 - 0.32e2 * t323 * t409 + 0.16e2 * t353 * t402 + 0.16e2 * t357 * t402) * t115 * t327) * (PREF_V_CFCA2_At + PREF_V_CA_At / 0.2e1) + (((0.64e2 * t172 * t175 * t182 * t23 * t286 * t84 - 0.128e3 * t172 * t19 * t362 * t365 + 0.64e2 * t172 * t175 * t226 + 0.64e2 * t198 * t199 * t422 + 0.64e2 * t198 * t199 * t426 + t421) * t230 + (-0.16e2 * t121 * t386 * t96 - 0.16e2 * t126 * t386 * t96 - 0.16e2 * t240 * t286 * t84 + 0.32e2 * t286 * t396 * t453 - 0.16e2 * t240 * t286 + 0.32e2 * t243 * t247) * t291 + (-0.16e2 * t294 * t322 * t84 + 0.32e2 * t322 * t409 * t453 - 0.16e2 * t294 * t322 - 0.32e2 * t296 * t300 + 0.16e2 * t402 * t422 + 0.16e2 * t402 * t426) * t327 - t421 - 0.128e3 * t337 * t422 - 0.128e3 * t337 * t426 - 0.64e2 * t381 * t132 + 0.32e2 * t333 * t85 + 0.32e2 * t333 * t27) * t94 + (-0.16e2 * FH0 * t182 * t286 * t39 * t514 * t88 - 0.16e2 * FH0 * t182 * t286 * t507 * t55 * t88 - 0.32e2 * t198 * t199 * t47 * t49 + 0.32e2 * t198 * t199 * t47 * t60) * t115 * t230 + (0.4e1 * FH0 * t180 * t248 * t546 + 0.4e1 * FH0 * t180 * t393 * t536 - 0.8e1 * t199 * t254 * t60 + 0.8e1 * t254 * t256) * t115 * t291 + (-0.4e1 * FH0 * t182 * t301 * t564 - 0.4e1 * FH0 * t182 * t378 * t559 + 0.8e1 * t305 * t307 * t60 - 0.8e1 * t305 * t308) * t115 * t327 + (0.32e2 * FH0 * (t329 + t500 - t170 + t331 - t178 + 0.4e1) * t26 * t378 + 0.32e2 * FH0 * (t329 - t500 - t170 + t331 + t178 + 0.4e1) * t26 * t301 + 0.64e2 * t193 * t160 - 0.64e2 * t193 * t163) * t115) * (PREF_V_CFCA2_Bt + PREF_V_CA_Bt / 0.2e1));
}
