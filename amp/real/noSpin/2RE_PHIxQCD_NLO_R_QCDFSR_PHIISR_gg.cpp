
#include "AMP_HEADER.h"

double Eval_R_FSR_ISR (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);



  double t1;
  double t10;
  double t1000;
  double t1004;
  double t1005;
  double t1006;
  double t101;
  double t1013;
  double t1016;
  double t1018;
  double t102;
  double t1020;
  double t1028;
  double t1029;
  double t103;
  double t1032;
  double t1033;
  double t1036;
  double t1039;
  double t104;
  double t105;
  double t1052;
  double t107;
  double t108;
  double t1080;
  double t1087;
  double t1089;
  double t109;
  double t1093;
  double t1099;
  double t110;
  double t1107;
  double t1109;
  double t1115;
  double t1123;
  double t113;
  double t1134;
  double t1137;
  double t1139;
  double t1142;
  double t1149;
  double t116;
  double t117;
  double t12;
  double t123;
  double t124;
  double t125;
  double t128;
  double t129;
  double t13;
  double t130;
  double t131;
  double t132;
  double t135;
  double t136;
  double t137;
  double t138;
  double t14;
  double t140;
  double t142;
  double t145;
  double t148;
  double t149;
  double t15;
  double t150;
  double t157;
  double t16;
  double t160;
  double t161;
  double t164;
  double t165;
  double t166;
  double t169;
  double t17;
  double t171;
  double t175;
  double t176;
  double t177;
  double t18;
  double t182;
  double t183;
  double t187;
  double t188;
  double t194;
  double t2;
  double t20;
  double t205;
  double t206;
  double t21;
  double t227;
  double t228;
  double t23;
  double t233;
  double t234;
  double t24;
  double t244;
  double t25;
  double t261;
  double t276;
  double t281;
  double t282;
  double t3;
  double t302;
  double t303;
  double t304;
  double t305;
  double t307;
  double t308;
  double t31;
  double t311;
  double t315;
  double t317;
  double t32;
  double t321;
  double t328;
  double t33;
  double t330;
  double t331;
  double t332;
  double t334;
  double t336;
  double t338;
  double t339;
  double t349;
  double t350;
  double t352;
  double t358;
  double t36;
  double t360;
  double t361;
  double t37;
  double t373;
  double t375;
  double t389;
  double t39;
  double t393;
  double t394;
  double t395;
  double t396;
  double t398;
  double t399;
  double t4;
  double t40;
  double t410;
  double t413;
  double t42;
  double t421;
  double t423;
  double t424;
  double t425;
  double t428;
  double t43;
  double t433;
  double t438;
  double t444;
  double t445;
  double t446;
  double t447;
  double t452;
  double t453;
  double t454;
  double t459;
  double t46;
  double t47;
  double t5;
  double t501;
  double t506;
  double t507;
  double t508;
  double t509;
  double t510;
  double t512;
  double t517;
  double t522;
  double t526;
  double t528;
  double t529;
  double t533;
  double t539;
  double t54;
  double t544;
  double t547;
  double t55;
  double t551;
  double t56;
  double t560;
  double t561;
  double t564;
  double t569;
  double t576;
  double t583;
  double t588;
  double t594;
  double t597;
  double t6;
  double t600;
  double t601;
  double t606;
  double t61;
  double t616;
  double t62;
  double t63;
  double t64;
  double t641;
  double t643;
  double t645;
  double t649;
  double t653;
  double t654;
  double t656;
  double t66;
  double t662;
  double t663;
  double t67;
  double t672;
  double t676;
  double t677;
  double t678;
  double t685;
  double t69;
  double t691;
  double t693;
  double t694;
  double t695;
  double t7;
  double t71;
  double t713;
  double t714;
  double t72;
  double t725;
  double t73;
  double t733;
  double t735;
  double t737;
  double t74;
  double t741;
  double t748;
  double t75;
  double t76;
  double t762;
  double t779;
  double t787;
  double t795;
  double t8;
  double t80;
  double t800;
  double t806;
  double t807;
  double t808;
  double t809;
  double t81;
  double t810;
  double t812;
  double t814;
  double t815;
  double t82;
  double t824;
  double t827;
  double t83;
  double t836;
  double t84;
  double t842;
  double t85;
  double t86;
  double t878;
  double t9;
  double t900;
  double t92;
  double t95;
  double t97;
  double t980;
  double t982;
  double t986;
  double t989;
  double t992;
  t1 = EPS_(k2, p1, p2, p3);
  t2 = sp(p2, p3);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = sp(p1, p3);
  t6 = 0.1e1 / t5;
  t7 = t1 * t6;
  t8 = t4 + t7;
  t9 = 0.96e2 * t8;
  t10 = sp(p1, p2);
  t12 = t4 * t5;
  t13 = 0.96e2 * t12;
  t14 = 0.64e2 * t1;
  t15 = t1 * t2;
  t16 = t15 * t6;
  t17 = 0.96e2 * t16;
  t18 = t1 * t5;
  t20 = 0.32e2 * t18 + 0.32e2 * t15;
  t21 = 0.1e1 / t10;
  t23 = t10 * t9 + t20 * t21 - t13 - t14 - t17;
  t24 = sp(k2, p2);
  t25 = 0.1e1 / t24;
  t31 = -t10 * t9 - t20 * t21 + t13 + t14 + t17;
  t32 = sp(k2, p1);
  t33 = 0.1e1 / t32;
  t36 = sp(k1, p3);
  t37 = 0.1e1 / t36;
  t39 = sp(k2, p3);
  t40 = 0.1e1 / t39;
  t42 = sp(k1, p2);
  t43 = 0.1e1 / t42;
  t46 = sp(k1, p1);
  t47 = 0.1e1 / t46;
  t54 = t2 * t2;
  t55 = t5 * t5;
  t56 = -t54 - t55;
  t61 = t10 * t10;
  t62 = t6 * t61;
  t63 = 0.128e3 * t62;
  t64 = 0.3e1 * t54;
  t66 = t2 + t5;
  t67 = 0.1e1 / t66;
  t69 = t6 * t10;
  t71 = 0.64e2 * (t64 + t55) * t67 * t69;
  t72 = t54 * t2;
  t73 = 0.3e1 * t72;
  t74 = t54 * t5;
  t75 = t55 * t2;
  t76 = t55 * t5;
  t80 = 0.64e2 * (t73 + t74 + t75 - t76) * t67 * t6;
  t81 = 0.64e2 * t54;
  t82 = 0.128e3 * t55;
  t83 = t2 * t5;
  t84 = 0.64e2 * t83;
  t85 = t6 * t72;
  t86 = 0.128e3 * t85;
  t92 = t3 + t6;
  t95 = t6 * t2;
  t97 = t3 * t5;
  t101 = t54 * t54;
  t102 = 0.3e1 * t101;
  t103 = t72 * t5;
  t104 = 0.16e2 * t103;
  t105 = t54 * t55;
  t107 = t76 * t2;
  t108 = 0.16e2 * t107;
  t109 = t55 * t55;
  t110 = 0.3e1 * t109;
  t113 = t67 * t3;
  t116 = t3 * t76;
  t117 = 0.128e3 * t116;
  t123 = t2 - t5;
  t124 = t123 * t3;
  t125 = t67 * t6;
  t128 = t124 * t125 * t61 * t10;
  t129 = 0.128e3 * t128;
  t130 = 0.5e1 * t74;
  t131 = 0.3e1 * t75;
  t132 = 0.3e1 * t76;
  t135 = t66 * t66;
  t136 = 0.1e1 / t135;
  t137 = t3 * t136;
  t138 = t137 * t61;
  t140 = 0.64e2 * (t73 - t130 - t131 - t132) * t6 * t138;
  t142 = 0.6e1 * t107;
  t145 = t137 * t10;
  t148 = 0.2e1 * t101;
  t149 = 0.3e1 * t107;
  t150 = 0.2e1 * t109;
  t157 = 0.128e3 * t61;
  t160 = t67 * t2 * t5 * t10;
  t161 = 0.256e3 * t160;
  t164 = t3 * t61;
  t165 = 0.128e3 * t164;
  t166 = 0.3e1 * t55;
  t169 = t3 * t10;
  t171 = 0.64e2 * (t54 + t166) * t67 * t169;
  t175 = 0.64e2 * (t72 - t74 - t75 - t132) * t67 * t3;
  t176 = 0.64e2 * t55;
  t177 = 0.128e3 * t54;
  t182 = 0.3e1 * t74;
  t183 = 0.5e1 * t75;
  t187 = 0.64e2 * (t73 + t182 + t183 - t132) * t6 * t138;
  t188 = 0.6e1 * t103;
  t194 = 0.3e1 * t103;
  t205 = 0.4e1 * t83;
  t206 = 0.5e1 * t55;
  t227 = 0.8e1 * t105;
  t228 = 0.10e2 * t107;
  t233 = 0.5e1 * t103;
  t234 = 0.7e1 * t107;
  t244 = 0.128e3 * (t64 + t205 + t166) * t67 * t10;
  t261 = 0.5e1 * t54;
  t276 = 0.10e2 * t103;
  t281 = 0.7e1 * t103;
  t282 = 0.5e1 * t107;
  t302 = 0.256e3 * t6;
  t303 = -t56;
  t304 = t303 * t67;
  t305 = t6 * t21;
  t307 = 0.256e3 * t304 * t305;
  t308 = t6 * t54;
  t311 = 0.1e1 / t61;
  t315 = t32 * t32;
  t317 = 0.256e3 * t3;
  t321 = t3 * t55;
  t328 = t125 * t10;
  t330 = 0.256e3 * t124 * t328;
  t331 = t123 * t303;
  t332 = t137 * t6;
  t334 = 0.256e3 * t331 * t332;
  t336 = t125 * t21;
  t338 = 0.128e3 * t331 * t3 * t336;
  t339 = 0.256e3 * t95;
  t349 = t303 * t3;
  t350 = t67 * t21;
  t352 = 0.256e3 * t349 * t350;
  t358 = t24 * t24;
  t360 = 0.256e3 * t97;
  t361 = t113 * t21;
  t373 = t303 * t303;
  t375 = t3 * t6;
  t389 = t315 * t32 * t2;
  t393 = -0.128e3 * t66;
  t394 = t393 * t311;
  t395 = 0.3e1 * t83;
  t396 = 0.2e1 * t55;
  t398 = (t54 - t395 - t396) * t67;
  t399 = t305 * t39;
  t410 = 0.4e1 * t54;
  t413 = t350 * t39;
  t421 = 0.2e1 * t54;
  t423 = (t421 - t395 - t55) * t3;
  t424 = t39 * t39;
  t425 = t125 * t424;
  t428 = 0.5e1 * t83;
  t433 = 0.2e1 * t74;
  t438 = t54 - t83 + t396;
  t444 = 0.64e2 * t62;
  t445 = 0.4e1 * t2;
  t446 = 0.4e1 * t5;
  t447 = 0.8e1 * t83;
  t452 = 0.4e1 * t55;
  t453 = 0.11e2 * t75;
  t454 = 0.5e1 * t76;
  t459 = 0.2e1 * t76;
  t501 = 0.2e1 * t83;
  t506 = 0.4e1 * t72;
  t507 = 0.4e1 * t74;
  t508 = 0.4e1 * t75;
  t509 = 0.4e1 * t76;
  t510 = 0.7e1 * t101;
  t512 = 0.7e1 * t109;
  t517 = 0.2e1 * t105;
  t522 = t124 * t67;
  t526 = 0.128e3 * t522 * t69 * t424 * t39;
  t528 = t124 * t125 * t61;
  t529 = 0.64e2 * t528;
  t533 = 0.128e3 * (t54 + t83 + t396) * t136 * t69;
  t539 = 0.64e2 * t128;
  t544 = t375 * t61;
  t547 = 0.16e2 * t105;
  t551 = t375 * t10;
  t560 = 0.64e2 * t61;
  t561 = 0.128e3 * t160;
  t564 = t358 * t24;
  t569 = (t421 + t395 - t55) * t3;
  t576 = (t54 + t395 - t396) * t3;
  t583 = 0.2e1 * t75;
  t588 = t421 - t83 + t55;
  t594 = 0.64e2 * t164;
  t597 = t67 * t10;
  t600 = 0.5e1 * t72;
  t601 = 0.11e2 * t74;
  t606 = 0.2e1 * t72;
  t616 = 0.128e3 * (t421 + t83 + t55) * t136 * t169;
  t641 = 0.128e3 * t25 * t311 * t389 + (t394 + (0.128e3 * t21 * t5 + 0.128e3 * t398 * t399 - 0.128e3 * t95 - 0.128e3) * t25) * t315 + (-t393 * t311 * t24 + 0.128e3 * (t410 + t83 - t55) * t3 * t413 - t339 - 0.384e3 - 0.128e3 * t97 + 0.128e3 * (t72 + t75 + t76) * t3 * t305 + (-0.128e3 * t423 * t425 + (0.64e2 * (t54 + t428 + t396) * t3 * t328 + 0.128e3 * t123 * (t72 + t433 + t131 + t76) * t332 - 0.64e2 * t2 * t438 * t336) * t39 + t444 + 0.32e2 * (-t445 - t446 + t64 + t447 + t55) * t67 * t69 - 0.32e2 * (-t410 - t452 + t72 + t182 + t453 + t454) * t67 * t6 - 0.32e2 * (t410 - t452 + t74 - t131 - t459) * t6 * t21) * t25) * t32 + t394 * t358 + (-0.128e3 * (t54 - t83 - t452) * t67 * t399 - 0.128e3 * t95 - t360 - 0.384e3 + 0.128e3 * (t72 + t74 + t76) * t3 * t305) * t24 + 0.128e3 * (t54 - 0.6e1 * t83 + t55) * t3 * t425 + (0.64e2 * (t261 + 0.14e2 * t83 + t206) * t3 * t328 - 0.128e3 * (t101 + t194 + t149 + t109) * t136 * t375 + 0.64e2 * t66 * t21) * t39 - 0.64e2 * t92 * t61 + 0.32e2 * (t445 + t446 + t261 + t501 + t206) * t3 * t69 - 0.32e2 * (t506 + t507 + t508 + t509 + t510 + t104 + 0.26e2 * t105 + t108 + t512) * t3 * t125 + 0.32e2 * (t506 - t507 - t508 + t509 + t148 + t194 - t517 + t149 + t150) * t3 * t305 + (t526 + (0.64e2 * t438 * t6 * t67 + t529 - t533) * t424 + (-t539 - 0.32e2 * (-t410 + t452 + t73 + 0.19e2 * t74 + 0.13e2 * t75 + t454) * t136 * t544 + 0.32e2 * (-t506 + t507 - t508 + t509 + t101 + t276 + t547 + 0.22e2 * t107 + t512) * t136 * t551 - 0.32e2 * (-t506 - t507 + t508 + t509 + t103 + t517 + t234 + t150) * t3 * t125) * t39 + t560 + t561) * t25 + (0.128e3 * t311 * t564 * t5 + (0.128e3 * t2 * t21 - 0.128e3 * t413 * t569 - 0.128e3 * t97 - 0.128e3) * t358 + (0.128e3 * t576 * t425 + (0.64e2 * (t421 + t428 + t55) * t3 * t328 - 0.128e3 * t123 * (t72 + t182 + t583 + t76) * t332 - 0.64e2 * t5 * t588 * t361) * t39 + t594 + 0.32e2 * (-t445 - t446 + t54 + t447 + t166) * t3 * t597 - 0.32e2 * (-t410 - t452 + t600 + t601 + t131 + t76) * t3 * t67 + 0.32e2 * (t410 - t452 + t606 + t182 - t75) * t3 * t21) * t24 - t526 + (0.64e2 * t3 * t588 * t67 - t529 - t616) * t424 + (t539 - 0.32e2 * (t410 - t452 + t600 + 0.13e2 * t74 + 0.19e2 * t75 + t132) * t136 * t544 + 0.32e2 * (t506 - t507 + t508 - t509 + t510 + 0.22e2 * t103 + t547 + t228 + t109) * t136 * t551 - 0.32e2 * (t506 + t507 - t508 - t509 + t148 + t281 + t517 + t107) * t3 * t125) * t39 + t560 + t561) * t33;
  t643 = t40 * t311;
  t645 = 0.128e3 * t643 * t389;
  t649 = -0.128e3 * t311 * t54 * t6 - t302 + t307;
  t653 = 0.128e3 * t398 * t305;
  t654 = 0.128e3 * t69;
  t656 = (t64 + t501 + t55) * t67;
  t662 = 0.384e3 * t83;
  t663 = -t86 - t662;
  t672 = 0.128e3 * t3 * t311 * t55 + t317 - t352;
  t676 = t123 * t67 * t21;
  t677 = 0.256e3 * t676;
  t678 = 0.128e3 * t169;
  t685 = -0.256e3 * t54 + 0.256e3 * t55;
  t691 = t125 * t39;
  t693 = 0.128e3 * t423 * t691;
  t694 = 0.128e3 * t528;
  t695 = 0.7e1 * t72;
  t713 = 0.15e2 * t74;
  t714 = 0.15e2 * t75;
  t725 = -0.256e3 * t72 - 0.256e3 * t75;
  t733 = 0.128e3 * t643 * t564 * t5;
  t735 = 0.128e3 * t569 * t350;
  t737 = -0.3e1 * t5 + t2;
  t741 = t662 + t117;
  t748 = 0.128e3 * t576 * t691;
  t762 = 0.7e1 * t55;
  t779 = 0.256e3 * t74 + 0.256e3 * t76;
  t787 = 0.128e3 * t522 * t69 * t424;
  t795 = 0.7e1 * t76;
  t800 = 0.4e1 * t105;
  t806 = 0.4e1 * t101;
  t807 = 0.8e1 * t103;
  t808 = 0.8e1 * t107;
  t809 = 0.4e1 * t109;
  t810 = t101 * t5;
  t812 = t72 * t55;
  t814 = t76 * t54;
  t815 = t109 * t2;
  t824 = 0.256e3 * t83 * t676;
  t827 = t54 + t83 + t55;
  t836 = -t645 + (t649 * t40 * t24 - t653 + (0.256e3 * t2 * t303 * t336 + t311 * t663 - 0.128e3 * t6 * t656 + t654) * t40) * t315 + (t672 * t40 * t358 + (t330 - t334 - t677 + (-0.64e2 * t21 * t303 * t6 + t311 * t685 + 0.128e3 * t349 * t67 - t678) * t40) * t24 + t693 - t694 + 0.64e2 * (t695 + t433 + t75 - t459) * t136 * t551 - 0.128e3 * t2 * (t64 + t83 + t452) * t136 * t6 + 0.192e3 * (t72 - t74 + t508 + t459) * t67 * t305 + (-t444 + 0.32e2 * (t445 + t446 + t64 + t205 + t206) * t67 * t69 - 0.32e2 * (t410 + t452 + t72 + t713 + t714 + t454) * t67 * t6 + 0.32e2 * (t506 + t507 - t508 - t509 + 0.11e2 * t103 + t517 + t282 - t150) * t67 * t305 + t725 * t311) * t40) * t32 + t733 + (-t735 + (0.64e2 * t303 * t361 * t737 + t311 * t741 + t360) * t40) * t358 + (t748 + 0.64e2 * (t2 - 0.5e1 * t5) * t3 * t597 + 0.64e2 * t123 * (t72 - t74 + t75 - t132) * t332 - 0.64e2 * (t600 + t601 - t508 + t459) * t3 * t350 + (t594 - 0.32e2 * (t445 + t446 + t54 + t205 + t762) * t3 * t597 - 0.32e2 * (-t410 - t452 + t73 - t601 - 0.7e1 * t75 - t454) * t3 * t67 + 0.32e2 * (t506 + t507 - t508 - t509 + 0.6e1 * t101 - t194 - 0.9e1 * t107 - t150) * t3 * t350 + t779 * t311) * t40) * t24 - t787 + (-t529 + t533 - 0.64e2 * (t2 - 0.2e1 * t5) * t737 * t125) * t39 + t539 - 0.32e2 * (t410 - t452 + t73 - t74 - t131 - t795) * t136 * t544 + 0.32e2 * (t506 - t507 + t508 - t509 + t101 + t188 - t800 - t142 - 0.5e1 * t109) * t136 * t551 - 0.32e2 * (-0.2e1 * t109 * t5 + t806 + t807 - t808 - t809 + 0.5e1 * t810 - 0.5e1 * t812 + t814 - 0.7e1 * t815) * t136 * t375 - t824 + (-t560 + 0.64e2 * t656 * t10 - 0.256e3 * t2 * t827 * t67 + 0.64e2 * (t410 - t452 + t73 + t74 - t459) * t21) * t40;
  t842 = 0.3e1 * t2 - t5;
  t878 = 0.7e1 * t54;
  t900 = t54 + t501 + t166;
  t980 = t645 + (-t649 * t40 * t24 + t653 + (-0.64e2 * t303 * t336 * t842 - t311 * t663 + t339) * t40) * t315 + (-t672 * t40 * t358 + (-t330 + t334 + t677 + (-0.64e2 * t21 * t349 + 0.128e3 * t304 * t6 - t311 * t685 - t654) * t40) * t24 - t693 - 0.64e2 * (0.5e1 * t2 - t5) * t67 * t69 + 0.64e2 * t123 * (t73 - t74 + t75 - t76) * t332 - 0.64e2 * (t606 - t507 + t453 + t454) * t67 * t305 + (t444 - 0.32e2 * (t445 + t446 + t878 + t205 + t55) * t67 * t69 + 0.32e2 * (t410 + t452 + t600 + 0.7e1 * t74 + t453 - t132) * t67 * t6 - 0.32e2 * (t506 + t507 - t508 - t509 + t148 + 0.9e1 * t103 + t149 - 0.6e1 * t109) * t67 * t305 - t725 * t311) * t40) * t32 - t733 + (t735 + (-0.128e3 * t3 * t67 * t900 + 0.256e3 * t303 * t361 * t5 - t311 * t741 + t678) * t40) * t358 + (-t748 + t694 - 0.64e2 * (t606 - t74 - t583 - t795) * t136 * t551 - 0.128e3 * t5 * (t410 + t83 + t166) * t137 + 0.192e3 * (t606 + t507 - t75 + t76) * t3 * t350 + (-t594 + 0.32e2 * (t445 + t446 + t261 + t205 + t166) * t3 * t597 - 0.32e2 * (t410 + t452 + t600 + t713 + t714 + t76) * t3 * t67 - 0.32e2 * (t506 + t507 - t508 - t509 + t148 - t233 - t517 - 0.11e2 * t107) * t3 * t350 - t779 * t311) * t40) * t24 + t787 + (t529 + t616 - 0.64e2 * t842 * (0.2e1 * t2 - t5) * t113) * t39 - t539 + 0.32e2 * (t410 - t452 + t695 + t182 + t75 - t132) * t136 * t544 - 0.32e2 * (t506 - t507 + t508 - t509 + 0.5e1 * t101 + t188 + t800 - t142 - t109) * t136 * t551 + 0.32e2 * (0.2e1 * t101 * t2 + t806 + t807 - t808 - t809 + 0.7e1 * t810 - t812 + 0.5e1 * t814 - 0.5e1 * t815) * t136 * t375 + t824 + (-t560 + 0.64e2 * t900 * t67 * t10 - 0.256e3 * t5 * t827 * t67 - 0.64e2 * (t410 - t452 + t606 - t75 - t132) * t21) * t40;
  t982 = (-t302 + t307 + (-0.128e3 * t5 - 0.128e3 * t308) * t311) * t40 * t315 + ((t302 + t317 + (-0.256e3 * t95 - 0.256e3 * t97) * t21 + (0.128e3 * t5 + 0.128e3 * t321 + 0.128e3 * t2 + 0.128e3 * t308) * t311) * t40 * t24 + t330 - t334 + t338 + (t339 - 0.256e3 - 0.256e3 * t331 * t336 + (-0.128e3 * t55 - 0.128e3 * t54 + 0.128e3 * t85 + 0.128e3 * t83) * t311) * t40) * t32 + (-t317 + t352 + (-0.128e3 * t321 - 0.128e3 * t2) * t311) * t40 * t358 + (-t330 + t334 - t338 + (-0.256e3 + t360 + 0.256e3 * t331 * t361 + (-0.128e3 * t54 - 0.128e3 * t55 + 0.128e3 * t83 + 0.128e3 * t116) * t311) * t40) * t24 - 0.256e3 * t349 * t328 + 0.256e3 * t373 * t136 * t375 - 0.128e3 * t373 * t3 * t336 + (0.128e3 * t21 * t303 + 0.256e3 * t10 - 0.256e3 * t304) * t40 + t641 * t37 + t836 * t43 + t980 * t47;
  t986 = -0.512e3 * t8;
  t989 = 0.256e3 * t12 + 0.256e3 * t16;
  t992 = (t21 * t986 + t311 * t989) * t40;
  t1000 = (-t21 * t986 - t311 * t989) * t40 * t24;
  t1004 = t4 * t55;
  t1005 = t1 * t54;
  t1006 = t1005 * t6;
  t1013 = t311 * t315;
  t1016 = t1 * t311;
  t1018 = t1 * t66;
  t1020 = t1018 * t3 * t399;
  t1028 = 0.256e3 * t4;
  t1029 = 0.256e3 * t7;
  t1032 = 0.256e3 * t4 * t6 * t424;
  t1033 = t1 * t123;
  t1036 = 0.128e3 * t1033 * t3 * t328;
  t1039 = 0.256e3 * t1033 * t827 * t332;
  t1052 = t125 * t3;
  t1080 = t3 * t21;
  t1087 = t1 * t40;
  t1089 = 0.256e3 * t1087 * t1013;
  t1093 = 0.256e3 * t1018 * t375 * t21;
  t1099 = 0.768e3 * t18 + 0.256e3 * t1006;
  t1107 = 0.256e3 * t1087 * t311 * t358;
  t1109 = 0.512e3 * t1;
  t1115 = 0.256e3 * t1004 + 0.768e3 * t15;
  t1123 = 0.256e3 * t4 * t6 * t39;
  t1134 = 0.64e2 * t8;
  t1137 = 0.128e3 * t1;
  t1139 = 0.256e3 * t1;
  t1142 = 0.9e1 * t83;
  t1149 = 0.512e3 * t1 * t55 + 0.512e3 * t1005;
  return(Bt_fH_re * ((t23 * t25 + t31 * t33) * t37 + t31 * t40 * t43 + t23 * t40 * t47) * PREF_R_CA + Bt_fA_re * ((0.256e3 * t21 * t56 + 0.512e3 * t2 + 0.512e3 * t5) * t40 + ((-t63 + t71 - t80 + (-t81 - t82 + t84 + t86) * t21) * t25 * t32 + 0.128e3 * t92 * t61 + (-0.192e3 * t95 - 0.384e3 - 0.192e3 * t97) * t10 + 0.64e2 * (t102 + t104 + 0.18e2 * t105 + t108 + t110) * t6 * t113 + (-t117 - 0.320e3 * t54 - t86 - 0.320e3 * t55 - 0.128e3 * t83) * t21 + ((t129 - t140 + 0.64e2 * (t102 - 0.2e1 * t103 - t142 - t110) * t6 * t145 - 0.64e2 * (t148 - t103 - t149 - t150) * t6 * t113) * t39 - t157 - t161) * t25 + ((-t165 + t171 + t175 + (-t176 + t117 + t84 - t177) * t21) * t24 + (-t129 + t187 - 0.64e2 * (t102 + t188 + 0.2e1 * t107 - t110) * t6 * t145 + 0.64e2 * (t148 + t194 + t107 - t150) * t6 * t113) * t39 - t157 - t161) * t33) * t37 + ((t63 - 0.64e2 * (t64 + t205 + t206) * t67 * t69 + 0.64e2 * (t73 + t130 + t75 + t132) * t67 * t6 + (-t84 + t82 - t86 - 0.192e3 * t54) * t21) * t40 * t32 + (-t165 + t171 + t175 + (t117 + t176 - 0.256e3 * t54 + t84) * t21) * t40 * t24 - t129 + t140 - 0.64e2 * (t102 - t188 - t227 - t228 - t110) * t6 * t145 + 0.64e2 * (t148 - t233 - t234 - t150) * t6 * t113 + 0.256e3 * t21 * t54 + (t157 - t244 + 0.512e3 * t54 + (-0.384e3 * t72 - 0.128e3 * t74 + 0.256e3 * t76) * t21) * t40) * t43 + ((-t63 + t71 - t80 + (t84 + t81 - 0.256e3 * t55 + t86) * t21) * t40 * t32 + (t165 - 0.64e2 * (t261 + t205 + t166) * t67 * t169 + 0.64e2 * (t73 + t74 + t183 + t132) * t67 * t3 + (-t117 + t177 - 0.192e3 * t55 - t84) * t21) * t40 * t24 + t129 - t187 + 0.64e2 * (t102 + t276 + t227 + t142 - t110) * t6 * t145 - 0.64e2 * (t148 + t281 + t282 - t150) * t6 * t113 + 0.256e3 * t21 * t55 + (t157 - t244 + 0.512e3 * t55 + (0.256e3 * t72 - 0.128e3 * t75 - 0.384e3 * t76) * t21) * t40) * t47) * PREF_R_CA + At_fH_re * t982 * PREF_R_CA + At_fA_re * (t992 * t32 + t1000 + ((-0.512e3 * t12 + 0.512e3 * t16) * t21 + (0.256e3 * t1004 - 0.256e3 * t18 + 0.256e3 * t15 - 0.256e3 * t1006) * t311) * t40 + (-0.256e3 * t1 * t25 * t1013 + (0.256e3 * t1016 + (0.256e3 * t1020 - 0.256e3 * t7) * t25) * t32 - 0.256e3 * t1016 * t24 + t1028 - t1029 + (-t1032 + (-0.128e3 * t1033 * t2 * t336 - t1036 + t1039) * t39 + 0.64e2 * t1 * (0.9e1 * t2 + t5) * t551 - 0.64e2 * t1 * (t410 + t447 + t452 + t695 + t130 + t183 - t76) * t1052 + 0.64e2 * t1 * t588 * t305) * t25 + (0.256e3 * t1016 * t358 + (-0.256e3 * t1020 + 0.256e3 * t4) * t24 + t1032 + (-0.128e3 * t1033 * t361 * t5 - t1036 + t1039) * t39 - 0.64e2 * t1 * (t2 + 0.9e1 * t5) * t551 - 0.64e2 * t1 * (-t410 - t447 - t452 + t72 - t130 - t183 - t795) * t1052 - 0.64e2 * t1 * t438 * t1080) * t33) * t37 + (t1089 + (t992 * t24 - t1093 + ((-0.512e3 * t1 - 0.512e3 * t16) * t21 + t1099 * t311) * t40) * t32 + t1107 + (-t1093 + (-t1028 - t1029 + (-0.384e3 * t12 - t1109 + 0.128e3 * t16) * t21 + t1115 * t311) * t40) * t24 + t1123 - t1036 + 0.256e3 * t1 * (t606 + t507 + t183 + t76) * t332 - 0.128e3 * t1 * (t73 + t130 + 0.6e1 * t75 + t459) * t6 * t361 + (t1134 * t10 + 0.64e2 * t12 - t1137 + t1028 + (-0.192e3 * t15 + t1139) * t6 + 0.64e2 * t1 * (t421 - t1142 - t762) * t305 + t1149 * t311) * t40) * t43 + (-t1089 + (t1000 + t1093 + (t1028 + t1029 + (-0.128e3 * t12 + t1109 + 0.384e3 * t16) * t21 - t1099 * t311) * t40) * t32 - t1107 + (t1093 + ((0.512e3 * t12 + 0.512e3 * t1) * t21 - t1115 * t311) * t40) * t24 - t1123 - t1036 - 0.256e3 * t1 * (t72 + t130 + t508 + t459) * t332 + 0.128e3 * t1 * (t606 + 0.6e1 * t74 + t183 + t132) * t6 * t361 + (-t1134 * t10 + 0.192e3 * t12 + t1137 - t1028 + (-0.64e2 * t15 - t1139) * t6 + 0.64e2 * t1 * (t878 + t1142 - t396) * t1080 - t1149 * t311) * t40) * t47) * PREF_R_CA);
}
