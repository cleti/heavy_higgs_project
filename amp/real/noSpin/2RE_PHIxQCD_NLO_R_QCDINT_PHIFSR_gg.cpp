
#include "AMP_HEADER.h"

double Eval_R_INT_FSR (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);



  double t1;
  double t10;
  double t102;
  double t11;
  double t111;
  double t117;
  double t14;
  double t140;
  double t141;
  double t142;
  double t143;
  double t144;
  double t145;
  double t146;
  double t147;
  double t151;
  double t152;
  double t154;
  double t157;
  double t158;
  double t16;
  double t163;
  double t165;
  double t17;
  double t170;
  double t172;
  double t174;
  double t175;
  double t177;
  double t178;
  double t180;
  double t185;
  double t187;
  double t188;
  double t189;
  double t19;
  double t190;
  double t193;
  double t198;
  double t2;
  double t20;
  double t208;
  double t209;
  double t210;
  double t214;
  double t221;
  double t223;
  double t227;
  double t230;
  double t231;
  double t232;
  double t234;
  double t242;
  double t243;
  double t246;
  double t25;
  double t252;
  double t254;
  double t259;
  double t26;
  double t263;
  double t265;
  double t277;
  double t278;
  double t279;
  double t28;
  double t281;
  double t285;
  double t289;
  double t291;
  double t292;
  double t293;
  double t295;
  double t304;
  double t309;
  double t311;
  double t313;
  double t317;
  double t32;
  double t322;
  double t325;
  double t329;
  double t33;
  double t34;
  double t340;
  double t342;
  double t349;
  double t350;
  double t351;
  double t359;
  double t362;
  double t37;
  double t378;
  double t381;
  double t384;
  double t387;
  double t39;
  double t391;
  double t4;
  double t408;
  double t411;
  double t416;
  double t422;
  double t425;
  double t428;
  double t433;
  double t441;
  double t443;
  double t446;
  double t448;
  double t455;
  double t456;
  double t46;
  double t464;
  double t467;
  double t468;
  double t469;
  double t475;
  double t477;
  double t48;
  double t483;
  double t49;
  double t499;
  double t5;
  double t50;
  double t505;
  double t510;
  double t519;
  double t52;
  double t53;
  double t534;
  double t536;
  double t538;
  double t544;
  double t582;
  double t583;
  double t59;
  double t598;
  double t60;
  double t609;
  double t611;
  double t615;
  double t618;
  double t62;
  double t630;
  double t631;
  double t633;
  double t635;
  double t637;
  double t64;
  double t640;
  double t642;
  double t643;
  double t644;
  double t645;
  double t647;
  double t65;
  double t651;
  double t653;
  double t68;
  double t680;
  double t684;
  double t69;
  double t690;
  double t692;
  double t693;
  double t695;
  double t697;
  double t7;
  double t70;
  double t700;
  double t702;
  double t705;
  double t708;
  double t71;
  double t710;
  double t713;
  double t716;
  double t73;
  double t754;
  double t756;
  double t758;
  double t761;
  double t765;
  double t770;
  double t772;
  double t774;
  double t783;
  double t787;
  double t790;
  double t792;
  double t798;
  double t8;
  double t837;
  double t840;
  double t841;
  double t842;
  double t843;
  double t844;
  double t851;
  double t852;
  double t862;
  double t869;
  double t870;
  double t873;
  double t882;
  double t883;
  double t885;
  double t895;
  double t898;
  double t9;
  double t90;
  double t91;
  double t910;
  double t92;
  double t926;
  double t93;
  double t94;
  double t941;
  t1 = sp(p2, p3);
  t2 = sp(p1, p3);
  t4 = 0.1e1 / (t1 + t2);
  t5 = sp(k1, k2);
  t7 = 0.1e1 / (0.1e1 + t5);
  t8 = t4 * t7;
  t9 = sp(k2, p3);
  t10 = 0.1e1 / t9;
  t11 = sp(k1, p3);
  t14 = 0.2e1 * t2;
  t16 = (t1 + t14) * t4;
  t17 = t7 * t10;
  t19 = t1 * t4;
  t20 = 0.1e1 / t11;
  t25 = sp(k1, p1);
  t26 = t25 * t25;
  t28 = t11 * t11;
  t32 = sp(k2, p1);
  t33 = t10 * t32;
  t34 = t8 * t33;
  t37 = 0.9e1 * t2;
  t39 = (t37 + t1) * t4;
  t46 = t8 * t9;
  t48 = t2 * t4;
  t49 = t48 * t7;
  t50 = 0.1280e4 * t49;
  t52 = 0.384e3 * t7 * t2;
  t53 = 0.256e3 * t7;
  t59 = t9 * t9;
  t60 = t8 * t59;
  t62 = 0.5e1 * t1;
  t64 = (t37 + t62) * t4;
  t65 = t7 * t9;
  t68 = 0.2e1 * t1;
  t69 = t2 * t2;
  t70 = 0.3e1 * t69;
  t71 = t1 * t2;
  t73 = (t14 + t68 + t70 - t71) * t4;
  t90 = 0.5e1 * t2;
  t91 = 0.3e1 * t1;
  t92 = t90 + t91;
  t93 = t2 * t92;
  t94 = t8 * t10;
  t102 = t32 * t32;
  t111 = 0.13e2 * t2;
  t117 = t2 * (t2 + t1 + t69);
  t140 = t48 * t10;
  t141 = -t140 + t4;
  t142 = 0.128e3 * t141;
  t143 = 0.1e1 / t32;
  t144 = t142 * t143;
  t145 = t9 * t4;
  t146 = t48 - t145;
  t147 = 0.128e3 * t146;
  t151 = sp(k1, p2);
  t152 = 0.1e1 / t151;
  t154 = t26 * t25;
  t157 = sp(k2, p2);
  t158 = 0.1e1 / t157;
  t163 = t4 * t10;
  t165 = -t141;
  t170 = 0.640e3 * t4;
  t172 = (t2 - t68) * t4;
  t174 = 0.128e3 * t172 * t10;
  t175 = 0.3e1 * t2;
  t177 = (-0.2e1 + t175 - t1) * t4;
  t178 = 0.2e1 * t69;
  t180 = (-t2 + t1 + t178 - t71) * t4;
  t185 = 0.512e3 * t145;
  t187 = (t68 + t90) * t4;
  t188 = 0.128e3 * t187;
  t189 = t59 * t4;
  t190 = 0.192e3 * t189;
  t193 = (-0.4e1 + 0.11e2 * t2 + t68) * t4;
  t198 = (-t175 - t1 + 0.4e1 * t69 + t71) * t4;
  t208 = t32 * t4;
  t209 = t158 * t10;
  t210 = t208 * t209;
  t214 = (-0.4e1 + t37) * t4;
  t221 = 0.384e3 * t4;
  t223 = 0.128e3 * t187 * t10;
  t227 = 0.384e3 * t145;
  t230 = 0.256e3 * (0.1e1 + t2 + t1) * t4;
  t231 = 0.5e1 * t69;
  t232 = 0.2e1 * t71;
  t234 = (-t175 - t1 + t231 + t232) * t4;
  t242 = t9 * t2;
  t243 = t242 * t4;
  t246 = (-t2 + t1 + t70) * t4;
  t252 = t143 * t10;
  t254 = -t252 * t48 + t163;
  t259 = (-t1 + t14) * t4;
  t263 = (-0.2e1 + t90) * t4;
  t265 = t2 * (-0.2e1 + t14 - t1);
  t277 = 0.128e3 * (-0.4e1 + t175 - t62) * t4;
  t278 = 0.6e1 * t1;
  t279 = 0.3e1 * t71;
  t281 = (t14 + t278 + t70 - t279) * t4;
  t285 = (t175 - t68) * t4;
  t289 = (-t2 + t91 + t178 - t279) * t4;
  t291 = t69 * t2;
  t292 = 0.2e1 * t291;
  t293 = t69 * t1;
  t295 = (-t14 - t68 - t69 + t279 + t292 - t293) * t4;
  t304 = (-0.4e1 + t111 + t68) * t4;
  t309 = 0.5e1 * t71;
  t311 = (-0.10e2 * t2 - t278 + 0.11e2 * t69 + t309) * t4;
  t313 = t59 * t9;
  t317 = (t14 - 0.1e1) * t4;
  t322 = (-t37 - t91 + 0.8e1 * t69 + t232) * t4;
  t325 = 0.9e1 * t69;
  t329 = (t14 + t68 - t325 - t309 + 0.5e1 * t291 + 0.2e1 * t293) * t4;
  t340 = 0.7e1 * t2;
  t342 = (-0.4e1 + t340) * t4;
  t349 = t102 * t4;
  t350 = t349 * t209;
  t351 = 0.512e3 * t350;
  t359 = (t175 + t1) * t4;
  t362 = (-t340 - t1 + t325 + t232) * t4;
  t378 = 0.384e3 * t189;
  t381 = (-0.4e1 + 0.15e2 * t2 - t68) * t4;
  t384 = 0.4e1 * t71;
  t387 = 0.128e3 * (t2 - t91 + t69 + t384) * t4;
  t391 = (t14 + t68 - t325 - t309 + 0.6e1 * t291 + 0.3e1 * t293) * t4;
  t408 = t59 * t2 * t4;
  t411 = (t1 + t70) * t4;
  t416 = (-t14 - t68 - t69 + t279 + 0.3e1 * t291) * t4;
  t422 = t208 * t10;
  t425 = t69 * t4;
  t428 = 0.128e3 * t48;
  t433 = t349 * t10;
  t441 = 0.128e3 * t243;
  t443 = t2 * (-0.4e1 + t340 - t68);
  t446 = 0.4e1 * t1;
  t448 = t2 * (t446 + t231);
  t455 = -t4 - t140;
  t456 = 0.128e3 * t455;
  t464 = 0.256e3 * t359 * t9;
  t467 = 0.128e3 * t408;
  t468 = t68 + t2;
  t469 = t2 * t468;
  t475 = 0.128e3 * t2 * (-t175 - t62 + t232) * t4;
  t477 = t2 * (-t14 - t68 + t232 + t291);
  t483 = t102 * t2 * t4;
  t499 = t2 * (-0.4e1 + t90);
  t505 = t2 * (-0.14e2 * t2 - t278 + 0.7e1 * t69 + t232);
  t510 = t2 * (t14 + t68 - 0.6e1 * t69 - t384 + t292 + t293);
  t519 = t48 * t209;
  t534 = t102 * t32;
  t536 = t534 * t4 * t209;
  t538 = 0.128e3 * t4;
  t544 = 0.128e3 * t145;
  t582 = t4 * t158;
  t583 = t534 * t2 * t582;
  t598 = t59 * t69 * t4;
  t609 = 0.1e1 / t25;
  t611 = (t143 * t147 * t20 + t144) * t152 * t154 + (-t142 * t158 - t147 * t158 * t20 + ((0.64e2 * t143 * t165 - 0.128e3 * t163) * t11 + t170 - t174 + (-0.128e3 * t10 * t180 - 0.128e3 * t145 + 0.128e3 * t177) * t143 + (-t185 + t188 + (-0.64e2 * t193 * t9 + t190 + 0.128e3 * t198) * t143) * t20) * t152) * t26 + ((-0.384e3 * t210 + (-0.64e2 * t10 * t214 + 0.576e3 * t4) * t158) * t11 + (t221 + t223) * t158 * t32 + (0.128e3 * t10 * t234 - t227 - t230) * t158 + (-0.128e3 * t172 * t158 * t32 + (-t190 + 0.576e3 * t243 - 0.128e3 * t246) * t158) * t20 + (0.64e2 * t254 * t28 + (-0.704e3 * t4 + 0.128e3 * t259 * t10 + (0.128e3 * t163 * t265 + 0.128e3 * t145 - 0.128e3 * t263) * t143) * t11 + (t221 - t174) * t32 + 0.192e3 * t145 + t277 - 0.128e3 * t281 * t10 + (-0.128e3 * t10 * t295 - 0.64e2 * t285 * t9 + 0.128e3 * t289) * t143 + ((-t227 + t188) * t32 + 0.448e3 * t189 - 0.128e3 * t304 * t9 + 0.128e3 * t311 + (-0.128e3 * t313 * t4 + 0.256e3 * t317 * t59 - 0.128e3 * t322 * t9 + 0.128e3 * t329) * t143) * t20) * t152) * t25 + (0.448e3 * t210 + (0.64e2 * t10 * t342 - t221) * t158) * t28 + (-t351 + (-0.128e3 * t10 * t304 + 0.192e3 * t4) * t158 * t32 + (-0.128e3 * t10 * t362 + 0.256e3 * t359) * t158) * t11 + (t170 + t223) * t158 * t102 + (0.128e3 * t10 * t311 - 0.704e3 * t145 + t277) * t158 * t32 + (0.128e3 * t10 * t391 - 0.64e2 * t381 * t9 + t378 - t387) * t158 + ((-0.128e3 * t145 - 0.128e3 * t172) * t158 * t102 + (0.128e3 * t259 * t9 + 0.64e2 * t189 - 0.128e3 * t281) * t158 * t32 + (0.256e3 * t411 * t9 - 0.384e3 * t408 - 0.128e3 * t416) * t158) * t20 + ((-0.192e3 * t422 + t221 - 0.384e3 * t140 + (-0.192e3 * t10 * t425 + t428) * t143) * t28 + (0.128e3 * t433 + (-t221 + 0.576e3 * t140) * t32 - 0.64e2 * t381 + 0.256e3 * t411 * t10 + (0.64e2 * t163 * t448 - 0.64e2 * t4 * t443 + t441) * t143) * t11 + t456 * t102 + (-0.128e3 * t10 * t246 + 0.576e3 * t145 - t230) * t32 - t378 + t464 - t387 - 0.128e3 * t416 * t10 + (0.64e2 * t145 * t469 - 0.128e3 * t163 * t477 - t467 - t475) * t143 + (0.128e3 * t483 + (-0.64e2 * t214 * t9 + 0.128e3 * t234) * t32 + 0.64e2 * t342 * t59 - 0.128e3 * t362 * t9 + 0.128e3 * t391 + (-0.128e3 * t2 * t313 * t4 - 0.64e2 * t145 * t505 + 0.64e2 * t189 * t499 + 0.128e3 * t4 * t510) * t143) * t20) * t152 + ((-0.128e3 * t210 - 0.128e3 * t519) * t28 * t11 + (0.192e3 * t350 + 0.256e3 * t317 * t209 * t32 + (0.64e2 * t163 * t499 - t428) * t158) * t28 + (-0.128e3 * t536 + (-0.64e2 * t10 * t193 - t538) * t158 * t102 + (-0.128e3 * t10 * t322 - 0.64e2 * t285 + t544) * t158 * t32 + (-0.64e2 * t163 * t505 + 0.64e2 * t4 * t469 + t441) * t158) * t11 - t456 * t158 * t534 + (0.128e3 * t10 * t198 - 0.64e2 * t145 + 0.128e3 * t177) * t158 * t102 + (0.128e3 * t10 * t329 - 0.128e3 * t263 * t9 + 0.128e3 * t289) * t158 * t32 + (-0.64e2 * t145 * t443 + 0.128e3 * t163 * t510 + t467 - t475) * t158 + (-0.128e3 * t583 + (0.64e2 * t243 - 0.128e3 * t180) * t158 * t102 + (0.128e3 * t145 * t265 - 0.128e3 * t295 - 0.64e2 * t408) * t158 * t32 + (0.64e2 * t145 * t448 - 0.128e3 * t4 * t477 - 0.192e3 * t598) * t158) * t20) * t609;
  t615 = 0.256e3 * t165;
  t618 = -0.256e3 * t146;
  t630 = 0.768e3 * t4;
  t631 = t468 * t4;
  t633 = 0.256e3 * t631 * t10;
  t635 = (-t1 + t2) * t4;
  t637 = (t2 + t1 + t69 - t71) * t4;
  t640 = 0.256e3 * t10 * t637 - 0.256e3 * t635;
  t642 = t68 + t175;
  t643 = t642 * t4;
  t644 = 0.256e3 * t643;
  t645 = 0.128e3 * t189;
  t647 = (t340 + t68) * t4;
  t651 = (-t2 - t1 + t70 + t71) * t4;
  t653 = 0.128e3 * t647 * t9 - t645 - 0.256e3 * t651;
  t680 = 0.256e3 * t359 * t10;
  t684 = 0.256e3 * t163 * t71 + 0.512e3 * t48;
  t690 = -0.256e3 * t10 * t643 - t630;
  t692 = 0.512e3 * t635;
  t693 = 0.1024e4 * t2;
  t695 = (0.512e3 - t693) * t10;
  t697 = t2 + t91;
  t700 = 0.256e3 * t2 * t697 * t4;
  t702 = t2 * (-t14 - t68 + t69 + t232);
  t705 = -0.256e3 * t163 * t702 + 0.384e3 * t243 + t700;
  t708 = 0.256e3 * t145 - 0.256e3 * t631;
  t710 = t2 * t642;
  t713 = t2 * (-t14 - t68 + t70 + t232);
  t716 = 0.256e3 * t145 * t710 - 0.256e3 * t4 * t713;
  t754 = 0.384e3 * t48;
  t756 = 0.256e3 * t710 * t163;
  t758 = t69 * t468;
  t761 = 0.128e3 * t163 * t758 + 0.384e3 * t425;
  t765 = 0.256e3 * t455;
  t770 = -0.256e3 * t10 * t651 + t544 - 0.256e3 * t635;
  t772 = 0.512e3 * t243;
  t774 = 0.256e3 * t713 * t163;
  t783 = 0.384e3 * t9 * t69 * t4 + 0.256e3 * t758 * t4 + (0.256e3 * t69 - 0.256e3 * t291) * t10;
  t787 = -t441 + 0.256e3 * t637;
  t790 = 0.256e3 * t71 * t145;
  t792 = 0.256e3 * t702 * t4;
  t798 = 0.128e3 * t145 * t758 - 0.256e3 * t291 + 0.128e3 * t598 + 0.256e3 * t69;
  t837 = (t143 * t20 * t618 + t143 * t615) * t152 * t154 + (t615 * t158 + t618 * t158 * t20 + ((0.256e3 * t163 + t144) * t11 - t630 - t633 + t640 * t143 + (t143 * t653 + t185 - t644) * t20) * t152) * t26 + ((t142 * t158 + 0.256e3 * t210) * t11 + (-t630 - t633) * t158 * t32 + t640 * t158 + ((t185 - t644) * t158 * t32 + t653 * t158) * t20 + (-0.128e3 * t254 * t28 + (t143 * t684 + 0.512e3 * t422 + t538 + t680) * t11 + t690 * t32 + t544 - t692 + t695 + t705 * t143 + (t143 * t716 + t32 * t708 + t464 - t645 - t693 + 0.512e3) * t20) * t152) * t25 + (-0.128e3 * t210 + 0.128e3 * t519) * t28 + (t351 + (t538 + t680) * t158 * t32 + t684 * t158) * t11 + t690 * t158 * t102 + (t544 - t692 + t695) * t158 * t32 + t705 * t158 + (t708 * t158 * t102 + (-t645 + t464 + 0.512e3 - t693) * t158 * t32 + t716 * t158) * t20 + ((0.128e3 * t252 * t425 - 0.128e3 * t422) * t28 + (t143 * t761 + 0.128e3 * t33 * t647 + 0.256e3 * t433 + t754 + t756) * t11 + t765 * t102 + t770 * t32 + t772 + t700 - t774 + t783 * t143 + (t143 * t798 + t32 * t787 + t467 + 0.256e3 * t483 + t790 - t792) * t20) * t152 + ((0.128e3 * t209 * t425 - 0.128e3 * t350) * t28 + (0.256e3 * t536 + 0.128e3 * t647 * t209 * t102 + (t754 + t756) * t158 * t32 + t761 * t158) * t11 + t765 * t158 * t534 + t770 * t158 * t102 + (t772 + t700 - t774) * t158 * t32 + t783 * t158 + (0.256e3 * t583 + t787 * t158 * t102 + (t467 + t790 - t792) * t158 * t32 + t798 * t158) * t20) * t609;
  t840 = EPS_(k1, k2, p1, p3);
  t841 = t840 * t4;
  t842 = 0.256e3 * t841;
  t843 = t697 * t4;
  t844 = t840 * t10;
  t851 = 0.256e3 * t9 * t840 * t4;
  t852 = t92 * t4;
  t862 = t582 * t10;
  t869 = t840 * t158;
  t870 = t869 * t20;
  t873 = t2 * t840;
  t882 = 0.256e3 * t19 * t840;
  t883 = 0.4e1 * t2;
  t885 = (t883 + t446 + t69 - t71) * t4;
  t895 = t242 * t841;
  t898 = (-t883 - t446 + t231 + t279) * t4;
  t910 = t32 * t840;
  t926 = t910 * t158;
  t941 = t840 * t143;
  return(At_fH_re * (((-0.256e3 * t10 * t11 * t8 + 0.256e3 * t19 * t20 * t7 + 0.256e3 * t16 * t17 + 0.256e3 * t8) * t26 + (0.192e3 * t8 * t10 * t28 + (-0.64e2 * t17 * t39 - 0.256e3 * t34 - 0.192e3 * t8) * t11 + 0.512e3 * t8 * t32 - 0.448e3 * t46 + t50 + (t52 - t53) * t10 + (-0.256e3 * t32 * t8 * t9 + 0.64e2 * t64 * t65 - 0.128e3 * t7 * t73 - 0.64e2 * t60) * t20) * t25 + (0.64e2 * t17 * t48 - 0.64e2 * t34) * t28 + ((0.64e2 * t17 * t64 - 0.448e3 * t8) * t32 - 0.832e3 * t49 + 0.64e2 * t93 * t94) * t11 + (0.256e3 * t17 * t19 + 0.256e3 * t8) * t102 + (-0.128e3 * t17 * t73 - 0.192e3 * t46 + t50) * t32 - 0.832e3 * t48 * t65 + 0.128e3 * t2 * (t91 + t111) * t8 - 0.256e3 * t117 * t94 + ((0.256e3 * t16 * t7 - 0.256e3 * t46) * t102 + (-0.64e2 * t39 * t65 + t52 - t53 + 0.192e3 * t60) * t32 + 0.64e2 * t48 * t7 * t59 + 0.64e2 * t93 * t46 - 0.256e3 * t117 * t8) * t20) * PREF_R_CA + t611 * PREF_R_CFCA2) + Bt_fA_re * t837 * PREF_R_CFCA2 + At_fA_re * (((-0.128e3 * t843 * t844 - t842) * t143 + (-0.128e3 * t840 * t852 + t851) * t143 * t20) * t152 * t25 - 0.256e3 * t11 * t840 * t862 + (0.128e3 * t359 * t844 + t842) * t158 - 0.128e3 * t635 * t870 + ((-0.256e3 * t163 * t873 + 0.256e3 * t841) * t143 * t11 - t842 + 0.128e3 * t635 * t844 + (0.128e3 * t844 * t885 + t882) * t143 + (t851 - 0.128e3 * t359 * t840 + (-0.256e3 * t4 * t59 * t840 - 0.128e3 * t840 * t898 + 0.512e3 * t895) * t143) * t20) * t152 + (0.256e3 * t28 * t840 * t862 + (-0.512e3 * t862 * t873 - 0.256e3 * t862 * t910) * t11 + (0.128e3 * t844 * t852 + t842) * t158 * t32 + (0.128e3 * t844 * t898 - t851 - t882) * t158 + (0.128e3 * t843 * t926 + (-0.128e3 * t840 * t885 + 0.256e3 * t895) * t158) * t20) * t609) * PREF_R_CFCA2 + Bt_fH_re * ((-0.192e3 * t10 * t941 - 0.192e3 * t20 * t941) * t152 * t25 + 0.192e3 * t869 * t10 + 0.192e3 * t870 + (-0.192e3 * t873 * t252 - 0.192e3 * t844 + (-0.192e3 * t143 * t873 - 0.192e3 * t840) * t20) * t152 + (0.192e3 * t873 * t209 + 0.192e3 * t910 * t209 + (0.192e3 * t158 * t873 + 0.192e3 * t926) * t20) * t609) * PREF_R_CFCA2);
}
