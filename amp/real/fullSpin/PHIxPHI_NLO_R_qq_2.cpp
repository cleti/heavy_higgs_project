
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_QQ (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t10;
  double t100;
  double t1003;
  double t1008;
  double t1011;
  double t1012;
  double t1017;
  double t1020;
  double t1025;
  double t1042;
  double t105;
  double t1054;
  double t106;
  double t107;
  double t108;
  double t1085;
  double t11;
  double t1118;
  double t112;
  double t1137;
  double t1163;
  double t1172;
  double t1179;
  double t119;
  double t120;
  double t1202;
  double t123;
  double t1235;
  double t124;
  double t1262;
  double t1269;
  double t127;
  double t1276;
  double t1279;
  double t1286;
  double t131;
  double t1322;
  double t1327;
  double t1363;
  double t138;
  double t139;
  double t1401;
  double t143;
  double t1430;
  double t146;
  double t148;
  double t149;
  double t15;
  double t152;
  double t153;
  double t154;
  double t159;
  double t16;
  double t160;
  double t161;
  double t164;
  double t165;
  double t168;
  double t169;
  double t17;
  double t170;
  double t178;
  double t179;
  double t182;
  double t183;
  double t186;
  double t187;
  double t190;
  double t191;
  double t194;
  double t195;
  double t198;
  double t199;
  double t2;
  double t202;
  double t205;
  double t208;
  double t21;
  double t213;
  double t216;
  double t217;
  double t22;
  double t221;
  double t23;
  double t230;
  double t231;
  double t234;
  double t235;
  double t238;
  double t239;
  double t240;
  double t241;
  double t245;
  double t248;
  double t249;
  double t250;
  double t256;
  double t257;
  double t260;
  double t261;
  double t266;
  double t27;
  double t273;
  double t278;
  double t28;
  double t281;
  double t284;
  double t287;
  double t29;
  double t290;
  double t293;
  double t3;
  double t30;
  double t300;
  double t307;
  double t312;
  double t317;
  double t318;
  double t319;
  double t322;
  double t33;
  double t331;
  double t332;
  double t333;
  double t335;
  double t337;
  double t339;
  double t34;
  double t341;
  double t342;
  double t344;
  double t347;
  double t349;
  double t35;
  double t352;
  double t354;
  double t363;
  double t372;
  double t377;
  double t380;
  double t385;
  double t39;
  double t390;
  double t4;
  double t40;
  double t404;
  double t41;
  double t411;
  double t42;
  double t425;
  double t428;
  double t43;
  double t432;
  double t438;
  double t444;
  double t447;
  double t448;
  double t451;
  double t454;
  double t46;
  double t461;
  double t47;
  double t471;
  double t474;
  double t475;
  double t478;
  double t48;
  double t49;
  double t5;
  double t50;
  double t501;
  double t505;
  double t51;
  double t514;
  double t518;
  double t521;
  double t525;
  double t529;
  double t54;
  double t543;
  double t546;
  double t549;
  double t55;
  double t550;
  double t553;
  double t556;
  double t560;
  double t561;
  double t564;
  double t568;
  double t578;
  double t58;
  double t587;
  double t59;
  double t598;
  double t599;
  double t6;
  double t60;
  double t604;
  double t609;
  double t61;
  double t610;
  double t634;
  double t635;
  double t638;
  double t64;
  double t643;
  double t653;
  double t656;
  double t659;
  double t663;
  double t667;
  double t67;
  double t686;
  double t699;
  double t7;
  double t70;
  double t71;
  double t721;
  double t724;
  double t725;
  double t728;
  double t729;
  double t732;
  double t735;
  double t738;
  double t74;
  double t741;
  double t744;
  double t748;
  double t75;
  double t752;
  double t755;
  double t756;
  double t757;
  double t762;
  double t763;
  double t766;
  double t769;
  double t770;
  double t773;
  double t774;
  double t78;
  double t783;
  double t79;
  double t801;
  double t804;
  double t805;
  double t810;
  double t813;
  double t814;
  double t82;
  double t823;
  double t83;
  double t832;
  double t849;
  double t852;
  double t853;
  double t856;
  double t859;
  double t86;
  double t860;
  double t863;
  double t866;
  double t867;
  double t87;
  double t870;
  double t873;
  double t874;
  double t877;
  double t88;
  double t882;
  double t884;
  double t909;
  double t91;
  double t920;
  double t94;
  double t942;
  double t95;
  double t963;
  double t98;
  double t984;
  double t989;
  double t99;
  double t990;
  double t993;
  double t994;
  t1 = sp(p1, p3);
  t2 = t1 * t1;
  t3 = At_Bt_fA2_De * t2;
  t4 = EPS_(k1, p3, s1, s2);
  t5 = sp(p1, p2);
  t6 = t5 * t5;
  t7 = t4 * t6;
  t10 = sp(p2, p3);
  t11 = t10 * t10;
  t15 = sp(s1, p3);
  t16 = At_Bt_fA2_De * t15;
  t17 = EPS_(k1, p1, p3, s2);
  t21 = sp(s2, p3);
  t22 = At_Bt_fA2_De * t21;
  t23 = EPS_(k1, p1, p3, s1);
  t27 = sp(k1, p3);
  t28 = At_Bt_fA2_De * t27;
  t29 = t6 * t5;
  t30 = t29 * t4;
  t33 = t27 * t27;
  t34 = At_Bt_fA2_De * t33;
  t35 = EPS_(k1, p1, s1, s2);
  t39 = sp(k1, p2);
  t40 = t28 * t39;
  t41 = EPS_(k1, k2, p3, s2);
  t42 = t41 * t15;
  t43 = t42 * t10;
  t46 = At_Bt_fA2_De * t5;
  t47 = t46 * t1;
  t48 = sp(k1, s2);
  t49 = t48 * t27;
  t50 = EPS_(k1, k2, p1, s1);
  t51 = t49 * t50;
  t54 = EPS_(k1, k2, p3, s1);
  t55 = t49 * t54;
  t58 = sp(s1, k2);
  t59 = t58 * t27;
  t60 = EPS_(k1, k2, p1, s2);
  t61 = t59 * t60;
  t64 = t59 * t41;
  t67 = t59 * t17;
  t70 = t27 * t50;
  t71 = t70 * t21;
  t74 = sp(s2, p1);
  t75 = t70 * t74;
  t78 = t27 * t54;
  t79 = t78 * t21;
  t82 = -0.16e2 * t11 * t7 * At_Bt_fA2_De + 0.16e2 * t16 * t17 * t6 - 0.16e2 * t22 * t23 * t6 + 0.16e2 * t34 * t35 * t6 - 0.24e2 * t28 * t30 - 0.16e2 * t3 * t7 + 0.32e2 * t40 * t43 - 0.8e1 * t47 * t51 + 0.24e2 * t47 * t55 + 0.8e1 * t47 * t61 - 0.24e2 * t47 * t64 + 0.16e2 * t47 * t67 - 0.24e2 * t47 * t71 + 0.8e1 * t47 * t75 + 0.56e2 * t47 * t79;
  t83 = t78 * t74;
  t86 = t27 * t60;
  t87 = sp(s1, p1);
  t88 = t86 * t87;
  t91 = t86 * t15;
  t94 = t27 * t41;
  t95 = t94 * t87;
  t98 = t46 * t39;
  t99 = t41 * t87;
  t100 = t99 * t10;
  t105 = At_Bt_fA2_De * t1;
  t106 = t105 * t48;
  t107 = t58 * t39;
  t108 = EPS_(k1, k2, p1, p3);
  t112 = t27 * t39;
  t119 = t39 * t54;
  t120 = t119 * t10;
  t123 = t39 * t108;
  t124 = t123 * t15;
  t127 = t39 * t23;
  t131 = t105 * t58;
  t138 = t39 * t41;
  t139 = t138 * t10;
  t143 = t39 * t10 * t17;
  t146 = -0.16e2 * t10 * t106 * t127 + 0.8e1 * t106 * t107 * t108 + 0.8e1 * t106 * t112 * t50 - 0.24e2 * t106 * t112 * t54 + 0.24e2 * t112 * t131 * t41 - 0.8e1 * t112 * t131 * t60 - 0.8e1 * t100 * t98 + 0.32e2 * t106 * t120 + 0.16e2 * t106 * t124 - 0.32e2 * t131 * t139 + 0.16e2 * t131 * t143 + 0.64e2 * t43 * t98 - 0.40e2 * t47 * t83 - 0.8e1 * t47 * t88 + 0.24e2 * t47 * t91 + 0.40e2 * t47 * t95;
  t148 = t48 * t39;
  t149 = t148 * t54;
  t152 = t46 * t27;
  t153 = t60 * t15;
  t154 = t153 * t10;
  t159 = t105 * t27;
  t160 = t39 * t50;
  t161 = t160 * t21;
  t164 = t54 * t74;
  t165 = t164 * t10;
  t168 = EPS_(k1, k2, s1, s2);
  t169 = At_Bt_fA2_De * t168;
  t170 = t169 * t5;
  t178 = t58 * t108;
  t179 = t178 * t10;
  t182 = t108 * t87;
  t183 = t182 * t10;
  t186 = t108 * t15;
  t187 = t186 * t10;
  t190 = t48 * t58;
  t191 = t190 * t108;
  t194 = t48 * t54;
  t195 = t194 * t10;
  t198 = t48 * t108;
  t199 = t198 * t87;
  t202 = t198 * t15;
  t205 = t46 * t48;
  t208 = -0.24e2 * t47 * t149 + 0.32e2 * t152 * t154 + 0.40e2 * t152 * t100 + 0.16e2 * t159 * t161 + 0.8e1 * t98 * t165 - 0.8e1 * t170 * t112 * t10 - 0.8e1 * t170 * t1 * t27 * t39 + 0.16e2 * t106 * t179 - 0.16e2 * t106 * t183 + 0.24e2 * t106 * t187 - 0.8e1 * t47 * t191 - 0.32e2 * t47 * t195 + 0.8e1 * t47 * t199 - 0.24e2 * t47 * t202 - 0.8e1 * t205 * t179;
  t213 = t105 * t39;
  t216 = At_Bt_fA2_De * t48;
  t217 = t216 * t58;
  t221 = t216 * t27;
  t230 = t1 * t5;
  t231 = t230 * t4;
  t234 = t10 * t5;
  t235 = t234 * t4;
  t238 = At_Bt_fA2_De * t39;
  t239 = t238 * t23;
  t240 = t21 * t1;
  t241 = t240 * t5;
  t245 = t21 * t10 * t5;
  t248 = At_Bt_fA2_De * t58;
  t249 = t248 * t27;
  t250 = t39 * t60;
  t256 = t50 * t21;
  t257 = t256 * t10;
  t260 = t54 * t21;
  t261 = t260 * t10;
  t266 = 0.8e1 * t10 * t123 * t217 + 0.8e1 * t10 * t160 * t221 - 0.8e1 * t10 * t249 * t250 + 0.16e2 * t187 * t216 * t39 - 0.16e2 * t120 * t221 + 0.16e2 * t139 * t249 + 0.16e2 * t165 * t40 + 0.8e1 * t183 * t205 - 0.32e2 * t187 * t205 - 0.56e2 * t213 * t43 - 0.16e2 * t231 * t40 - 0.16e2 * t235 * t40 + 0.32e2 * t239 * t241 + 0.32e2 * t239 * t245 + 0.16e2 * t257 * t40 - 0.32e2 * t261 * t40;
  t273 = t234 * t35;
  t278 = t107 * t41;
  t281 = t107 * t17;
  t284 = t119 * t21;
  t287 = t119 * t74;
  t290 = t138 * t87;
  t293 = t138 * t15;
  t300 = t46 * t58;
  t307 = -0.16e2 * t100 * t40 - 0.24e2 * t120 * t205 - 0.16e2 * t124 * t205 + 0.16e2 * t139 * t300 - 0.8e1 * t143 * t300 - 0.16e2 * t152 * t161 - 0.16e2 * t154 * t40 - 0.56e2 * t159 * t235 + 0.16e2 * t159 * t273 + 0.16e2 * t278 * t47 - 0.8e1 * t281 * t47 - 0.56e2 * t284 * t47 + 0.8e1 * t287 * t47 - 0.8e1 * t290 * t47 + 0.56e2 * t293 * t47;
  t312 = t250 * t15;
  t317 = t2 * t1;
  t318 = At_Bt_fA2_De * t317;
  t319 = t5 * t4;
  t322 = t11 * t10;
  t331 = t317 * t5;
  t332 = sp(s1, s2);
  t333 = t332 * At2_fH2_De;
  t335 = t2 * t6;
  t337 = t322 * t5;
  t339 = t11 * t6;
  t341 = t2 * t10;
  t342 = t5 * At2_fH2_De;
  t344 = t2 * t5;
  t347 = t1 * t11;
  t349 = t11 * t5;
  t352 = 0.8e1 * t10 * t30 * At_Bt_fA2_De + 0.8e1 * t319 * t322 * At_Bt_fA2_De + 0.8e1 * t105 * t30 + 0.32e2 * t152 * t284 - 0.16e2 * t152 * t287 + 0.16e2 * t152 * t312 + 0.8e1 * t318 * t319 - t331 * t333 + t333 * t335 - t333 * t337 + t333 * t339 - 0.2e1 * t333 * t344 - 0.2e1 * t333 * t349 - 0.16e2 * t34 * t7 + t341 * t342 + t342 * t347;
  t354 = t332 * At2_fA2_De;
  t363 = t5 * At2_fA2_De;
  t372 = t332 * Bt2_fH2_De;
  t377 = t5 * Bt2_fH2_De;
  t380 = t332 * Bt2_fA2_De;
  t385 = -0.4e1 * t331 * t354 + t331 * t372 + 0.4e1 * t331 * t380 + 0.4e1 * t335 * t354 - t335 * t372 - 0.4e1 * t335 * t380 - 0.4e1 * t337 * t354 + t337 * t372 + 0.4e1 * t339 * t354 - t339 * t372 + 0.4e1 * t341 * t363 + t341 * t377 - 0.8e1 * t344 * t354 + 0.4e1 * t347 * t363 + t347 * t377 - 0.8e1 * t349 * t354;
  t390 = t5 * Bt2_fA2_De;
  t404 = t168 * t5;
  t411 = t58 * t41;
  t425 = t23 * t21;
  t428 = -0.2e1 * t11 * t404 * At_Bt_fH2_De + 0.8e1 * t16 * t41 * t6 + 0.8e1 * t169 * t33 * t6 + 0.8e1 * t17 * t318 * t58 - 0.2e1 * t2 * t404 * At_Bt_fH2_De - 0.8e1 * t22 * t54 * t6 - 0.8e1 * t164 * t318 + 0.8e1 * t260 * t318 - 0.8e1 * t318 * t411 - 0.8e1 * t318 * t42 - 0.8e1 * t318 * t425 + 0.8e1 * t318 * t99 + 0.4e1 * t337 * t380 - 0.4e1 * t339 * t380 + 0.4e1 * t341 * t390 + 0.4e1 * t347 * t390;
  t432 = t23 * t74;
  t438 = t15 * t17;
  t444 = t322 * t17;
  t447 = At_Bt_fA2_De * t54;
  t448 = t21 * t322;
  t451 = t74 * t322;
  t454 = At_Bt_fA2_De * t41;
  t461 = At_Bt_fA2_De * t23;
  t471 = t234 * t48;
  t474 = t461 * t27;
  t475 = t230 * t48;
  t478 = 0.48e2 * t1 * t461 * t471 - 0.16e2 * t15 * t322 * t454 - 0.8e1 * t17 * t318 * t87 - 0.8e1 * t248 * t322 * t41 + 0.8e1 * t322 * t454 * t87 + 0.16e2 * t152 * t290 - 0.32e2 * t152 * t293 + 0.8e1 * t248 * t444 + 0.8e1 * t318 * t432 + 0.8e1 * t318 * t438 + 0.16e2 * t447 * t448 - 0.8e1 * t447 * t451 - 0.16e2 * t448 * t461 + 0.8e1 * t451 * t461 - 0.32e2 * t474 * t475;
  t501 = t411 * t10;
  t505 = t58 * t10 * t17;
  t514 = 0.16e2 * t100 * t213 - 0.32e2 * t100 * t47 - 0.32e2 * t159 * t284 + 0.16e2 * t159 * t287 - 0.16e2 * t159 * t290 + 0.32e2 * t159 * t293 - 0.16e2 * t159 * t312 - 0.16e2 * t165 * t213 + 0.32e2 * t165 * t47 + 0.16e2 * t183 * t221 - 0.16e2 * t187 * t221 + 0.56e2 * t213 * t261 - 0.80e2 * t261 * t47 + 0.32e2 * t47 * t501 - 0.32e2 * t47 * t505 - 0.24e2 * t471 * t474;
  t518 = t425 * t10;
  t521 = t432 * t10;
  t525 = t87 * t10 * t17;
  t529 = t15 * t10 * t17;
  t543 = t70 * t10;
  t546 = t78 * t10;
  t549 = t27 * t108;
  t550 = t549 * t87;
  t553 = t549 * t15;
  t556 = t27 * t23;
  t560 = 0.24e2 * t10 * t106 * t556 - 0.16e2 * t106 * t108 * t59 + 0.16e2 * t106 * t543 - 0.40e2 * t106 * t546 + 0.16e2 * t106 * t550 - 0.16e2 * t106 * t553 - 0.32e2 * t152 * t518 + 0.32e2 * t152 * t521 - 0.32e2 * t152 * t525 + 0.32e2 * t152 * t529 + 0.80e2 * t43 * t47 + 0.80e2 * t47 * t518 - 0.32e2 * t47 * t521 + 0.32e2 * t47 * t525 - 0.80e2 * t47 * t529;
  t561 = t86 * t10;
  t564 = t94 * t10;
  t568 = t27 * t10 * t17;
  t578 = t48 * t23;
  t587 = At_Bt_fA2_De * t29;
  t598 = t50 * t74;
  t599 = t598 * t10;
  t604 = -0.8e1 * t216 * t23 * t322 + 0.8e1 * t216 * t322 * t54 - 0.8e1 * t444 * t87 * At_Bt_fA2_De - 0.16e2 * t131 * t561 + 0.40e2 * t131 * t564 - 0.24e2 * t131 * t568 + 0.24e2 * t159 * t257 - 0.56e2 * t159 * t261 - 0.16e2 * t159 * t599 + 0.16e2 * t16 * t444 + 0.8e1 * t194 * t318 - 0.16e2 * t260 * t587 - 0.8e1 * t318 * t578 + 0.16e2 * t42 * t587 + 0.16e2 * t425 * t587 - 0.16e2 * t438 * t587;
  t609 = t60 * t87;
  t610 = t609 * t10;
  t634 = t238 * t15;
  t635 = t230 * t17;
  t638 = t234 * t17;
  t643 = -0.8e1 * t10 * t217 * t549 - 0.48e2 * t100 * t159 - 0.64e2 * t152 * t43 - 0.24e2 * t154 * t159 + 0.48e2 * t159 * t165 + 0.56e2 * t159 * t43 + 0.32e2 * t159 * t518 - 0.32e2 * t159 * t521 + 0.32e2 * t159 * t525 - 0.32e2 * t159 * t529 + 0.16e2 * t159 * t610 + 0.8e1 * t239 * t475 - 0.64e2 * t261 * t98 - 0.32e2 * t634 * t635 - 0.32e2 * t634 * t638;
  t653 = t94 * t15;
  t656 = t556 * t21;
  t659 = t556 * t74;
  t663 = t27 * t87 * t17;
  t667 = t27 * t15 * t17;
  t686 = 0.8e1 * t239 * t471 - 0.32e2 * t239 * t240 * t10 + 0.32e2 * t634 * t1 * t10 * t17 - 0.56e2 * t47 * t653 - 0.32e2 * t47 * t656 + 0.32e2 * t47 * t659 - 0.32e2 * t47 * t663 + 0.32e2 * t47 * t667 - 0.8e1 * t205 * t543 + 0.8e1 * t205 * t546 - 0.16e2 * t205 * t550 + 0.16e2 * t205 * t553 + 0.8e1 * t300 * t561 - 0.16e2 * t300 * t564 + 0.8e1 * t300 * t568 - 0.32e2 * t152 * t257;
  t699 = t1 * t6;
  t721 = t411 * t11;
  t724 = t58 * t11;
  t725 = t724 * t17;
  t728 = -0.24e2 * t28 * t344 * t4 + 0.8e1 * t28 * t349 * t35 - 0.8e1 * t28 * t35 * t699 + 0.48e2 * t28 * t4 * t699 + 0.24e2 * t100 * t3 - 0.24e2 * t105 * t721 + 0.24e2 * t105 * t725 - 0.40e2 * t152 * t165 + 0.64e2 * t152 * t261 + 0.8e1 * t152 * t599 - 0.8e1 * t152 * t610 - 0.32e2 * t3 * t43 - 0.32e2 * t3 * t518 + 0.24e2 * t3 * t521 - 0.24e2 * t3 * t525 + 0.32e2 * t3 * t529;
  t729 = t260 * t11;
  t732 = t164 * t11;
  t735 = t99 * t11;
  t738 = t42 * t11;
  t741 = t425 * t11;
  t744 = t432 * t11;
  t748 = t87 * t11 * t17;
  t752 = t15 * t11 * t17;
  t755 = At_Bt_fA2_De * t6;
  t756 = t1 * t48;
  t757 = t756 * t54;
  t762 = t1 * t58;
  t763 = t762 * t41;
  t766 = t762 * t17;
  t769 = t1 * t54;
  t770 = t769 * t74;
  t773 = t10 * t6;
  t774 = t773 * t4;
  t783 = 0.24e2 * t344 * t461 * t48 - 0.16e2 * t461 * t48 * t699 + 0.40e2 * t105 * t729 - 0.24e2 * t105 * t732 + 0.24e2 * t105 * t735 - 0.40e2 * t105 * t738 - 0.40e2 * t105 * t741 + 0.24e2 * t105 * t744 - 0.24e2 * t105 * t748 + 0.40e2 * t105 * t752 - 0.32e2 * t105 * t774 + 0.8e1 * t195 * t755 - 0.8e1 * t46 * t763 - 0.8e1 * t46 * t766 - 0.8e1 * t46 * t770 + 0.8e1 * t755 * t757;
  t801 = t108 * t21 * t87;
  t804 = t108 * t74;
  t805 = t804 * t15;
  t810 = t769 * t21;
  t813 = t1 * t41;
  t814 = t813 * t87;
  t823 = At_Bt_fA2_De * t108;
  t832 = 0.24e2 * t461 * t349 * t48 - 0.16e2 * t461 * t773 * t48 + 0.16e2 * t46 * t261 - 0.8e1 * t46 * t165 - 0.16e2 * t46 * t43 + 0.16e2 * t46 * t801 - 0.16e2 * t46 * t805 - 0.8e1 * t105 * t199 + 0.8e1 * t46 * t810 + 0.8e1 * t46 * t814 + 0.8e1 * t46 * t100 - 0.16e2 * t105 * t801 + 0.16e2 * t105 * t805 - 0.16e2 * t823 * t21 * t87 * t10 + 0.16e2 * t823 * t74 * t15 * t10;
  t849 = t578 * t10;
  t852 = t33 * t50;
  t853 = t852 * t21;
  t856 = t852 * t74;
  t859 = t33 * t54;
  t860 = t859 * t21;
  t863 = t859 * t74;
  t866 = t33 * t60;
  t867 = t866 * t87;
  t870 = t866 * t15;
  t873 = t33 * t41;
  t874 = t873 * t87;
  t877 = t873 * t15;
  t882 = 0.8e1 * t10 * t248 * t804 + 0.8e1 * t105 * t178 * t21 + 0.8e1 * t105 * t178 * t74 + 0.8e1 * t23 * t46 * t756 + 0.8e1 * t195 * t46 + 0.24e2 * t3 * t505 + 0.8e1 * t46 * t757 + 0.8e1 * t46 * t849 + 0.16e2 * t46 * t853 - 0.16e2 * t46 * t856 - 0.16e2 * t46 * t860 + 0.16e2 * t46 * t863 + 0.16e2 * t46 * t867 - 0.16e2 * t46 * t870 - 0.16e2 * t46 * t874 + 0.16e2 * t46 * t877;
  t884 = t349 * t4;
  t909 = t27 * t6;
  t920 = -0.16e2 * t230 * t34 * t35 + 0.8e1 * t28 * t344 * t35 - 0.8e1 * t28 * t35 * t773 + 0.16e2 * t461 * t48 * t909 + 0.24e2 * t105 * t884 + 0.16e2 * t231 * t34 + 0.24e2 * t235 * t3 + 0.16e2 * t235 * t34 - 0.16e2 * t238 * t741 + 0.16e2 * t238 * t752 + 0.16e2 * t241 * t461 + 0.16e2 * t245 * t461 - 0.16e2 * t273 * t34 + 0.56e2 * t28 * t774 - 0.32e2 * t28 * t884;
  t942 = t813 * t15;
  t963 = -0.16e2 * t16 * t635 - 0.16e2 * t16 * t638 + 0.16e2 * t238 * t15 * t6 * t17 - 0.16e2 * t238 * t21 * t6 * t23 + 0.16e2 * t238 * t909 * t4 - 0.8e1 * t105 * t202 - 0.8e1 * t216 * t183 - t344 * t190 * At2_fH2_De - 0.8e1 * t46 * t942 - 0.8e1 * t46 * t501 - 0.8e1 * t46 * t505 - 0.8e1 * t755 * t1 * t87 * t17 + 0.40e2 * t755 * t1 * t15 * t17 - 0.8e1 * t755 * t501 + 0.8e1 * t755 * t505 + 0.48e2 * t755 * t261;
  t984 = t39 * t39;
  t989 = t984 * t54;
  t990 = t989 * t21;
  t993 = t984 * t41;
  t994 = t993 * t15;
  t1003 = At_Bt_fA2_De * t984;
  t1008 = -0.8e1 * t755 * t165 + 0.8e1 * t755 * t100 - 0.48e2 * t755 * t43 - 0.48e2 * t755 * t518 + 0.8e1 * t755 * t521 - 0.8e1 * t755 * t525 + 0.48e2 * t755 * t529 - 0.8e1 * t169 * t5 * t33 * t10 - 0.8e1 * t105 * t58 * t984 * t41 + 0.16e2 * t105 * t990 - 0.16e2 * t105 * t994 + 0.8e1 * t216 * t989 * t10 - 0.8e1 * t248 * t993 * t10 + 0.16e2 * t1003 * t261 - 0.16e2 * t1003 * t43;
  t1011 = t5 * t332;
  t1012 = t1011 * Bt2_fA2_De;
  t1017 = t1011 * At2_fH2_De;
  t1020 = t1011 * At2_fA2_De;
  t1025 = t1011 * Bt2_fH2_De;
  t1042 = 0.4e1 * t1012 * t341 + 0.4e1 * t1012 * t347 - t1017 * t341 - t1017 * t347 - 0.4e1 * t1020 * t341 - 0.4e1 * t1020 * t347 + t1025 * t341 + t1025 * t347 - 0.16e2 * t46 * t735 + 0.48e2 * t46 * t738 + 0.48e2 * t46 * t741 - 0.16e2 * t46 * t744 + 0.16e2 * t46 * t748 - 0.48e2 * t46 * t752 + 0.8e1 * t55 * t755 - 0.8e1 * t755 * t763;
  t1054 = t1 * t23;
  t1085 = 0.8e1 * t755 * t766 + 0.40e2 * t755 * t810 - 0.8e1 * t755 * t770 + 0.8e1 * t755 * t814 - 0.40e2 * t755 * t942 - 0.40e2 * t755 * t1054 * t21 + 0.8e1 * t755 * t1054 * t74 - 0.8e1 * t248 * t27 * t11 * t17 + 0.16e2 * t28 * t256 * t11 - 0.8e1 * t28 * t598 * t11 - 0.32e2 * t28 * t729 + 0.24e2 * t28 * t732 + 0.8e1 * t105 * t48 * t984 * t54 - 0.24e2 * t3 * t501 + 0.16e2 * t46 * t721 - 0.16e2 * t46 * t725;
  t1118 = 0.24e2 * t3 * t64 + 0.24e2 * t3 * t653 + 0.16e2 * t3 * t656 - 0.16e2 * t3 * t659 + 0.16e2 * t3 * t663 - 0.16e2 * t3 * t667 - 0.16e2 * t3 * t67 + 0.8e1 * t3 * t71 - 0.8e1 * t3 * t75 - 0.24e2 * t3 * t79 + 0.24e2 * t3 * t83 + 0.8e1 * t3 * t88 - 0.8e1 * t3 * t91 - 0.24e2 * t3 * t95 - 0.48e2 * t46 * t729 + 0.16e2 * t46 * t732;
  t1137 = t5 * t48;
  t1163 = -0.8e1 * t10 * t216 * t852 + 0.8e1 * t10 * t216 * t859 + 0.8e1 * t10 * t248 * t866 - 0.8e1 * t10 * t248 * t873 + 0.8e1 * t11 * t216 * t556 + 0.8e1 * t11 * t216 * t70 - 0.16e2 * t11 * t216 * t78 - 0.8e1 * t11 * t248 * t86 + 0.16e2 * t11 * t248 * t94 - 0.4e1 * t1137 * t724 * At2_fA2_De + 0.16e2 * t105 * t874 - 0.16e2 * t105 * t877 - 0.16e2 * t257 * t34 + 0.16e2 * t261 * t34 + 0.16e2 * t34 * t599;
  t1172 = t48 * t33;
  t1179 = t58 * t33;
  t1202 = -0.16e2 * t105 * t1172 * t50 + 0.16e2 * t105 * t1172 * t54 - 0.16e2 * t105 * t1179 * t41 + 0.16e2 * t105 * t1179 * t60 + 0.16e2 * t100 * t34 - 0.16e2 * t105 * t853 + 0.16e2 * t105 * t856 + 0.16e2 * t105 * t860 - 0.16e2 * t105 * t863 - 0.16e2 * t105 * t867 + 0.16e2 * t105 * t870 + 0.16e2 * t154 * t34 - 0.16e2 * t165 * t34 - 0.16e2 * t34 * t610 + 0.16e2 * t71 * t755 - 0.32e2 * t755 * t79;
  t1235 = 0.16e2 * t23 * t3 * t49 + 0.32e2 * t238 * t729 - 0.8e1 * t238 * t732 + 0.8e1 * t238 * t735 + 0.8e1 * t3 * t51 - 0.24e2 * t3 * t55 - 0.8e1 * t3 * t61 + 0.32e2 * t653 * t755 + 0.16e2 * t656 * t755 - 0.16e2 * t659 * t755 + 0.16e2 * t663 * t755 - 0.16e2 * t667 * t755 + 0.16e2 * t755 * t83 - 0.16e2 * t755 * t91 - 0.16e2 * t755 * t95;
  t1262 = t2 * t58;
  t1269 = t2 * t54;
  t1276 = -0.16e2 * t11 * t153 * t28 + 0.8e1 * t11 * t28 * t609 - 0.16e2 * t1262 * t17 * t46 + 0.16e2 * t1262 * t41 * t46 - 0.32e2 * t1269 * t21 * t46 + 0.16e2 * t1269 * t46 * t74 - 0.32e2 * t238 * t738 - 0.24e2 * t28 * t735 + 0.32e2 * t28 * t738 + 0.16e2 * t28 * t741 - 0.16e2 * t28 * t744 + 0.16e2 * t28 * t748 - 0.16e2 * t28 * t752 - 0.16e2 * t34 * t43 - 0.16e2 * t46 * t990 + 0.16e2 * t46 * t994;
  t1279 = t2 * t41;
  t1286 = t2 * t23;
  t1322 = -0.16e2 * t46 * t1279 * t87 + 0.32e2 * t46 * t1279 * t15 + 0.32e2 * t46 * t1286 * t21 - 0.16e2 * t46 * t1286 * t74 + 0.16e2 * t46 * t2 * t87 * t17 - 0.32e2 * t46 * t2 * t15 * t17 + t344 * t190 * Bt2_fH2_De + t1137 * t724 * Bt2_fH2_De + 0.4e1 * t344 * t190 * Bt2_fA2_De + 0.4e1 * t1137 * t724 * Bt2_fA2_De - t1137 * t724 * At2_fH2_De - 0.4e1 * t344 * t190 * At2_fA2_De + 0.32e2 * t3 * t261 - 0.24e2 * t3 * t165 + 0.16e2 * t755 * t202;
  t1327 = t194 * t11;
  t1363 = -0.16e2 * t46 * t2 * t48 * t54 - 0.16e2 * t46 * t1327 + 0.8e1 * t3 * t191 + 0.24e2 * t3 * t195 - 0.8e1 * t3 * t199 + 0.8e1 * t3 * t202 - 0.24e2 * t3 * t849 + 0.24e2 * t105 * t1327 - 0.24e2 * t105 * t578 * t11 + 0.8e1 * t216 * t178 * t11 - 0.8e1 * t216 * t182 * t11 + 0.16e2 * t216 * t186 * t11 + 0.32e2 * t755 * t284 - 0.32e2 * t755 * t293 + 0.16e2 * t3 * t149 - 0.8e1 * t3 * t148 * t23;
  t1401 = -0.16e2 * t3 * t278 + 0.8e1 * t3 * t281 + 0.24e2 * t3 * t284 - 0.8e1 * t3 * t287 + 0.8e1 * t3 * t290 - 0.24e2 * t3 * t293 - 0.16e2 * t3 * t127 * t21 + 0.16e2 * t3 * t39 * t15 * t17 + 0.16e2 * t216 * t119 * t11 - 0.8e1 * t216 * t127 * t11 - 0.16e2 * t248 * t138 * t11 + 0.8e1 * t248 * t39 * t11 * t17 + t331 * At2_fH2_De - t335 * At2_fH2_De + t337 * At2_fH2_De - t339 * At2_fH2_De;
  t1430 = 0.4e1 * t331 * At2_fA2_De + 0.4e1 * t331 * Bt2_fA2_De + t331 * Bt2_fH2_De - 0.4e1 * t335 * At2_fA2_De - 0.4e1 * t335 * Bt2_fA2_De - t335 * Bt2_fH2_De + 0.4e1 * t337 * At2_fA2_De + 0.4e1 * t337 * Bt2_fA2_De + t337 * Bt2_fH2_De - 0.4e1 * t339 * At2_fA2_De - 0.4e1 * t339 * Bt2_fA2_De - t339 * Bt2_fH2_De + 0.8e1 * t344 * At2_fA2_De + 0.2e1 * t344 * At2_fH2_De + 0.8e1 * t349 * At2_fA2_De + 0.2e1 * t349 * At2_fH2_De;
  return(-0.128e3 * PREF_V_CA * (t428 + t478 + t514 + t560 + t604 + t643 + t882 + t920 + t963 + t82 + t208 + t686 + t728 + t783 + t832 + t266 + t307 + t352 + t385 + t1363 + t1401 + t1276 + t1322 + t146 + t1430 + t1008 + t1042 + t1235 + t1085 + t1202 + t1118 + t1163) / t6);
}
