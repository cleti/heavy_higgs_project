

double Eval_R_FSR_FSR (AMP_2_3_ARGS)
{

  AMP_2_3_4VEC_REFS(ps);


  double t1;
  double t10;
  double t100;
  double t101;
  double t106;
  double t107;
  double t11;
  double t112;
  double t113;
  double t115;
  double t119;
  double t12;
  double t120;
  double t121;
  double t122;
  double t124;
  double t128;
  double t132;
  double t134;
  double t135;
  double t136;
  double t137;
  double t139;
  double t14;
  double t141;
  double t143;
  double t146;
  double t147;
  double t148;
  double t149;
  double t15;
  double t157;
  double t158;
  double t16;
  double t160;
  double t162;
  double t163;
  double t164;
  double t166;
  double t169;
  double t17;
  double t170;
  double t171;
  double t172;
  double t174;
  double t176;
  double t185;
  double t187;
  double t188;
  double t189;
  double t190;
  double t194;
  double t195;
  double t196;
  double t2;
  double t200;
  double t201;
  double t202;
  double t206;
  double t208;
  double t213;
  double t219;
  double t224;
  double t228;
  double t23;
  double t234;
  double t238;
  double t239;
  double t24;
  double t240;
  double t241;
  double t243;
  double t245;
  double t247;
  double t250;
  double t251;
  double t252;
  double t254;
  double t257;
  double t258;
  double t259;
  double t26;
  double t260;
  double t266;
  double t270;
  double t272;
  double t273;
  double t277;
  double t279;
  double t283;
  double t286;
  double t292;
  double t293;
  double t299;
  double t3;
  double t30;
  double t301;
  double t308;
  double t309;
  double t31;
  double t310;
  double t317;
  double t319;
  double t327;
  double t336;
  double t338;
  double t34;
  double t355;
  double t359;
  double t360;
  double t361;
  double t364;
  double t367;
  double t370;
  double t387;
  double t388;
  double t39;
  double t390;
  double t399;
  double t400;
  double t412;
  double t413;
  double t43;
  double t44;
  double t45;
  double t450;
  double t454;
  double t456;
  double t460;
  double t461;
  double t47;
  double t472;
  double t473;
  double t475;
  double t48;
  double t481;
  double t486;
  double t490;
  double t491;
  double t492;
  double t494;
  double t496;
  double t498;
  double t5;
  double t501;
  double t503;
  double t505;
  double t507;
  double t51;
  double t516;
  double t518;
  double t520;
  double t521;
  double t525;
  double t526;
  double t528;
  double t532;
  double t543;
  double t548;
  double t55;
  double t550;
  double t558;
  double t560;
  double t566;
  double t568;
  double t575;
  double t576;
  double t579;
  double t585;
  double t586;
  double t587;
  double t589;
  double t59;
  double t598;
  double t6;
  double t60;
  double t600;
  double t603;
  double t608;
  double t61;
  double t610;
  double t613;
  double t615;
  double t619;
  double t621;
  double t622;
  double t624;
  double t627;
  double t632;
  double t634;
  double t636;
  double t64;
  double t643;
  double t646;
  double t659;
  double t660;
  double t662;
  double t665;
  double t666;
  double t673;
  double t675;
  double t679;
  double t681;
  double t685;
  double t690;
  double t7;
  double t711;
  double t718;
  double t72;
  double t743;
  double t749;
  double t75;
  double t750;
  double t753;
  double t762;
  double t764;
  double t77;
  double t792;
  double t796;
  double t799;
  double t8;
  double t801;
  double t806;
  double t810;
  double t811;
  double t815;
  double t824;
  double t826;
  double t83;
  double t843;
  double t881;
  double t9;
  double t90;
  double t919;
  double t95;
  double t97;
  double t99;
  t1 = EPS_(k1, k2, p1, p3);
  t2 = sp(p2, p3);
  t3 = sp(p1, p3);
  t5 = 0.1e1 / (t2 + t3);
  t6 = t1 * t5;
  t7 = sp(k2, p3);
  t8 = 0.1e1 / t7;
  t9 = sp(k2, p2);
  t10 = 0.1e1 / t9;
  t11 = t8 * t10;
  t12 = t6 * t11;
  t14 = sp(k2, p1);
  t15 = 0.1e1 / t14;
  t16 = t8 * t15;
  t17 = t6 * t16;
  t23 = sp(k1, p3);
  t24 = 0.1e1 / t23;
  t26 = t6 * t8;
  t30 = sp(k1, p2);
  t31 = 0.1e1 / t30;
  t34 = sp(k1, p1);
  t39 = t11 * t14;
  t43 = 0.256e3 * t6;
  t44 = 0.5e1 * t3;
  t45 = 0.3e1 * t2;
  t47 = t1 * (t44 + t45);
  t48 = t5 * t8;
  t51 = 0.128e3 * t47 * t48 - t43;
  t55 = t10 * t14;
  t59 = t6 * t7;
  t60 = t2 - t3;
  t61 = t1 * t60;
  t64 = 0.128e3 * t5 * t61 + 0.128e3 * t59;
  t72 = 0.128e3 * t6 * t8 * t23;
  t75 = 0.256e3 * t6 * t8 * t14;
  t77 = 0.128e3 * t61 * t48;
  t83 = -0.256e3 * t14 * t6 - 0.128e3 * t47 * t5 + 0.384e3 * t59;
  t90 = 0.1e1 / t34;
  t95 = t48 * t15;
  t97 = t48 * t10;
  t99 = t8 + t5;
  t100 = 0.256e3 * t99;
  t101 = t10 * t100;
  t106 = t8 - t5;
  t107 = 0.256e3 * t106;
  t112 = t34 * t34;
  t113 = t112 * t34;
  t115 = 0.128e3 * t95;
  t119 = t48 * t55;
  t120 = 0.768e3 * t119;
  t121 = 0.512e3 * t48;
  t122 = 0.11e2 * t3;
  t124 = (t45 + t122) * t5;
  t128 = (t44 + t2) * t5;
  t132 = 0.768e3 * t99 * t10;
  t134 = 0.512e3 * t5;
  t135 = 0.1024e4 * t8;
  t136 = t5 * t7;
  t137 = 0.384e3 * t136;
  t139 = (t45 - t44) * t5;
  t141 = 0.768e3 * t3;
  t143 = (-0.256e3 + t141) * t8;
  t146 = 0.128e3 * t136;
  t147 = 0.3e1 * t3;
  t148 = t2 - t147;
  t149 = t148 * t5;
  t157 = t23 * t23;
  t158 = 0.1e1 / t157;
  t160 = t48 * t23;
  t162 = t48 * t14;
  t163 = 0.768e3 * t162;
  t164 = 0.7e1 * t3;
  t166 = (-t164 + t2) * t5;
  t169 = t7 * t7;
  t170 = 0.1e1 / t169;
  t171 = 0.256e3 * t170;
  t172 = 0.768e3 * t106;
  t174 = 0.9e1 * t3;
  t176 = (t2 + t174) * t5;
  t185 = t48 * t10 * t157;
  t187 = 0.128e3 * t48;
  t188 = 0.10e2 * t3;
  t189 = t45 + t188;
  t190 = t189 * t5;
  t194 = (0.128e3 * t190 * t8 - 0.128e3 * t5) * t10;
  t195 = t15 * t3;
  t196 = t48 * t195;
  t200 = t14 * t14;
  t201 = t10 * t200;
  t202 = t48 * t201;
  t206 = 0.7e1 * t2;
  t208 = (t206 + t122) * t5;
  t213 = t3 * (t206 + 0.17e2 * t3);
  t219 = t3 * (t2 + t147);
  t224 = 0.768e3 * t136;
  t228 = (-0.512e3 + 0.1536e4 * t3) * t8;
  t234 = (t147 + t206) * t5;
  t238 = (-0.768e3 + 0.2304e4 * t3) * t8;
  t239 = t5 * t169;
  t240 = 0.128e3 * t239;
  t241 = 0.6e1 * t3;
  t243 = (t2 - t241) * t5;
  t245 = 0.128e3 * t243 * t7;
  t247 = t3 * (-t147 + t206);
  t250 = 0.512e3 * t3;
  t251 = t3 * t3;
  t252 = 0.768e3 * t251;
  t254 = (-t250 + t252) * t8;
  t257 = t136 * t3;
  t258 = 0.128e3 * t257;
  t259 = t3 * t60;
  t260 = t259 * t5;
  t266 = t10 * t3;
  t270 = t48 * t157;
  t272 = 0.128e3 * t5;
  t273 = t243 * t8;
  t277 = t48 * t200;
  t279 = t149 * t8;
  t283 = t259 * t48;
  t286 = 0.512e3 * t170 * t3;
  t292 = 0.128e3 * t190 * t7;
  t293 = t219 * t5;
  t299 = -0.128e3 * t185 + (t120 + t187 + t194 + 0.128e3 * t196) * t23 - 0.768e3 * t202 - 0.512e3 * t128 * t39 - 0.128e3 * t208 * t8 - t171 + (-0.128e3 * t213 * t48 + t146 + 0.512e3) * t10 - 0.128e3 * t219 * t95 + (t132 * t200 + (0.1536e4 * t8 + (-t224 - 0.512e3 * t149 + t228) * t10) * t14 - t146 - 0.128e3 * t234 + t238 + (-0.128e3 * t247 * t5 + t240 + t245 + t254) * t10 + (-t258 - 0.128e3 * t260 + t254) * t15) * t24 + (-0.512e3 * t195 - 0.512e3 * t55 - 0.512e3 * t266 - 0.512e3) * t158 + (0.128e3 * t270 + (-t163 + t272 + 0.128e3 * t273) * t23 + 0.768e3 * t277 + (-0.512e3 * t279 - 0.512e3 * t170) * t14 - t146 + 0.512e3 - 0.640e3 * t283 - t286 + (t172 * t200 + (t224 - 0.512e3 * t128 + t228) * t14 - t240 + t292 - 0.640e3 * t293 + t254) * t24) * t31;
  t301 = t48 * t266;
  t308 = t5 * t3;
  t309 = t164 + t45;
  t310 = t3 * t309;
  t317 = t200 * t14;
  t319 = t48 * t10 * t317;
  t327 = 0.512e3 * t170;
  t336 = 0.2e1 * t3;
  t338 = t251 * (t336 + t2);
  t355 = 0.1024e4 * t3;
  t359 = t239 * t3;
  t360 = 0.128e3 * t359;
  t361 = t3 * t148;
  t364 = t251 * t2;
  t367 = t251 * t3;
  t370 = (-0.256e3 * t251 + 0.256e3 * t367) * t8;
  t387 = t3 * t48 + t162;
  t388 = 0.128e3 * t387;
  t390 = 0.384e3 * t277;
  t399 = t48 * t317;
  t400 = 0.256e3 * t399;
  t412 = 0.256e3 * t170 * t251;
  t413 = t107 * t317;
  t450 = (-0.256e3 * t95 - 0.256e3 * t97 + (t100 * t15 + t101) * t24 + (t107 * t24 + 0.256e3 * t48) * t31) * t113 + ((t115 + 0.384e3 * t97) * t23 - t120 - t121 - 0.128e3 * t124 * t11 - 0.128e3 * t128 * t16 + (t132 * t14 + t134 + t135 + (-t137 - 0.128e3 * t139 + t143) * t10 + (-t146 - 0.128e3 * t149 + t143) * t15) * t24 + (-0.256e3 * t10 - 0.256e3 * t15) * t158 + (-0.384e3 * t160 + t163 - 0.128e3 * t166 * t8 - t171 + (t14 * t172 + t137 + t143 - 0.128e3 * t176) * t24) * t31) * t112 + t299 * t34 + (-0.128e3 * t301 - 0.128e3 * t119) * t157 + (0.384e3 * t202 + (-t187 + t194) * t14 + (0.128e3 * t310 * t48 - 0.128e3 * t308) * t10) * t23 - 0.256e3 * t319 + (-0.128e3 * t11 * t176 + t121) * t200 + (-0.128e3 * t234 * t8 - t327 + (-0.640e3 * t219 * t48 + t146 + 0.512e3) * t10) * t14 - t286 - 0.768e3 * t8 * t3 + (-0.512e3 * t338 * t48 + t250 + t258) * t10 + (t101 * t317 + (t135 - t134 + (-t137 - 0.128e3 * t166 + t143) * t10) * t200 + (t146 - 0.128e3 * t208 + t238 + (t240 + t245 - 0.640e3 * t260 + t254) * t10) * t14 - t141 + (-t355 + 0.1536e4 * t251) * t8 + (0.128e3 * t136 * t361 - 0.512e3 * t364 * t5 + t360 + t370) * t10 + t370 * t15) * t24 + (-0.256e3 * t201 + (-0.256e3 - 0.512e3 * t266) * t14 - t250 - 0.256e3 * t10 * t251 - 0.256e3 * t15 * t251) * t158 + (t388 * t157 + (-t390 + (0.128e3 * t5 + 0.128e3 * t273) * t14 + 0.128e3 * t308 + 0.128e3 * t361 * t48) * t23 + t400 + (-0.128e3 * t139 * t8 - t171) * t200 + (-0.128e3 * t247 * t48 - t146 - t286 + 0.512e3) * t14 - t258 + t250 - 0.512e3 * t364 * t48 - t412 + (t413 + (t137 - 0.128e3 * t124 + t143) * t200 + (-0.128e3 * t213 * t5 - t240 + t254 + t292) * t14 - t360 + 0.128e3 * t310 * t136 - 0.512e3 * t338 * t5 + t370) * t24) * t31 + ((-0.128e3 * t14 * t3 * t48 - 0.128e3 * t277) * t23 + t400 + (-0.128e3 * t279 - t171) * t200 + (-0.128e3 * t283 - t286) * t14 - t412 + (t413 + (t146 - 0.128e3 * t128 + t143) * t200 + (t258 - 0.128e3 * t293 + t254) * t14 + t370) * t24) * t90;
  t454 = 0.128e3 * t99;
  t456 = -t454;
  t460 = -t106;
  t461 = 0.128e3 * t460;
  t472 = 0.768e3 * t48;
  t473 = 0.640e3 * t5;
  t475 = (0.4e1 + t45 + t3) * t5;
  t481 = (-0.4e1 + t2 + t164) * t5;
  t486 = -t99;
  t490 = 0.512e3 * t8;
  t491 = 0.768e3 * t5;
  t492 = 0.64e2 * t136;
  t494 = (0.4e1 + t2 + t164) * t5;
  t496 = 0.1e1 + t3;
  t498 = 0.128e3 * t496 * t8;
  t501 = 0.192e3 * t136;
  t503 = (0.4e1 + t2 - t44) * t5;
  t505 = 0.1e1 - t3;
  t507 = 0.384e3 * t505 * t8;
  t516 = 0.640e3 * t162;
  t518 = (0.4e1 + t2 - t174) * t5;
  t520 = 0.64e2 * t518 * t8;
  t521 = 0.128e3 * t170;
  t525 = 0.5e1 * t2;
  t526 = 0.15e2 * t3;
  t528 = (-0.4e1 + t525 + t526) * t5;
  t532 = (0.384e3 - 0.640e3 * t3) * t8;
  t543 = (-0.4e1 + t2 + t44) * t5;
  t548 = 0.192e3 * t5;
  t550 = (-0.2e1 + t147) * t5;
  t558 = 0.13e2 * t3;
  t560 = (-0.4e1 + t45 + t558) * t5;
  t566 = 0.1152e4 * t5;
  t568 = (-0.4e1 + t206 + t558) * t5;
  t575 = 0.64e2 * (0.29e2 * t3 + t525) * t5;
  t576 = t3 * t2;
  t579 = (t2 - t147 + t576 + 0.6e1 * t251) * t5;
  t585 = 0.64e2 * t309 * t5;
  t586 = 0.2e1 * t576;
  t587 = 0.5e1 * t251;
  t589 = (-t2 - t44 + t586 + t587) * t5;
  t598 = 0.640e3 * t136;
  t600 = (0.4e1 + t45 - t147) * t5;
  t603 = 0.1024e4 * t505 * t8;
  t608 = 0.576e3 * t136;
  t610 = (0.4e1 + t206 - t526) * t5;
  t613 = 0.1152e4 * t505 * t8;
  t615 = (0.2e1 + t2 + t3) * t5;
  t619 = (t2 + t44 + t586 + t587) * t5;
  t621 = 0.256e3 * t3;
  t622 = 0.128e3 * t251;
  t624 = (-0.256e3 + t621 - t622) * t8;
  t627 = 0.64e2 * t239;
  t632 = (-t2 + t147 + t576 - 0.2e1 * t251) * t5;
  t634 = 0.384e3 * t251;
  t636 = (-0.256e3 + t141 - t634) * t8;
  t643 = -0.256e3 * t496;
  t646 = -0.256e3 * t505;
  t659 = 0.64e2 * (t558 + t525) * t5;
  t660 = 0.3e1 * t576;
  t662 = (-t45 + t3 + t660) * t5;
  t665 = -0.256e3 + t250;
  t666 = t665 * t170;
  t673 = 0.8e1 * t3;
  t675 = (-0.2e1 + t2 + t673) * t5;
  t679 = 0.11e2 * t251;
  t681 = (-t45 - t164 + 0.4e1 * t576 + t679) * t5;
  t685 = (-0.256e3 + t355 - 0.640e3 * t251) * t8;
  t690 = -0.64e2 * t185 + (-0.640e3 * t119 - 0.704e3 * t48 + (-0.64e2 * t543 * t8 + t134) * t10 + (-0.128e3 * t550 * t8 + t548) * t15) * t23 + 0.640e3 * t202 + (0.128e3 * t560 * t8 - t134) * t10 * t14 + t566 + 0.64e2 * t568 * t8 - t521 + (0.128e3 * t579 * t8 + 0.448e3 * t136 - t575) * t10 + (0.128e3 * t589 * t8 + t146 - t585) * t15 + (0.640e3 * t486 * t10 * t200 + (-0.768e3 * t8 + (t598 + 0.128e3 * t600 + t603) * t10) * t14 + t608 + 0.64e2 * t610 + t613 + (-0.128e3 * t615 * t7 - t240 + 0.128e3 * t619 + t624) * t10 + (-0.64e2 * t503 * t7 - t627 + 0.128e3 * t632 + t636) * t15) * t24 + (t10 * t643 + t15 * t646 + 0.512e3 * t55 + 0.512e3) * t158 + (-0.64e2 * t270 + (t516 - t272 - t520) * t23 - t390 + (0.128e3 * t600 * t8 - t134 + t327) * t14 + t501 - t659 + 0.128e3 * t662 * t8 + t666 + (0.384e3 * t460 * t200 + (-t598 + 0.128e3 * t560 + t603) * t14 + 0.256e3 * t239 - 0.128e3 * t675 * t7 + 0.128e3 * t681 + t685) * t24) * t31;
  t711 = t3 * (-0.4e1 + 0.4e1 * t2 + t122);
  t718 = t3 * (-0.4e1 + t45 + t673);
  t743 = (-0.10e2 * t2 - t188 + t660 - t251) * t5;
  t749 = 0.128e3 * t3 * t189 * t5;
  t750 = 0.2e1 * t2;
  t753 = t3 * (t750 - t336 + 0.5e1 * t576 + t679);
  t762 = 0.128e3 * t3 * (t147 + t750) * t5;
  t764 = t3 * (-t750 - t241 + t660 + t587);
  t792 = t3 * (0.4e1 + t750 + t3);
  t796 = t3 * (t750 + t241 + t660 + t587);
  t799 = 0.128e3 * t367;
  t801 = (-t621 + t622 - t799) * t8;
  t806 = t3 * (0.4e1 + t2 - t336);
  t810 = t3 * (-0.2e1 + t3);
  t811 = t60 * t5;
  t815 = (-t621 + t634 - t799) * t8;
  t824 = -t621 - t622;
  t826 = -t621 + t622;
  t843 = 0.128e3 * t399;
  t881 = 0.64e2 * t503 * t8;
  t919 = (-0.128e3 * t97 + t115 + (t10 * t454 + t15 * t456) * t24 + (t24 * t461 - t187) * t31) * t113 + ((-0.64e2 * t95 + 0.192e3 * t97) * t23 + 0.384e3 * t119 + t472 + (-0.64e2 * t475 * t8 - t473) * t10 + (0.64e2 * t481 * t8 - t272) * t15 + (0.384e3 * t486 * t10 * t14 - t490 - t491 + (-t492 + 0.64e2 * t494 + t498) * t10 + (t501 + 0.64e2 * t503 + t507) * t15) * t24 + (0.128e3 * t15 - 0.128e3 * t10) * t158 + (0.192e3 * t160 - t516 + t272 + t520 + t521 + (0.640e3 * t14 * t460 - 0.320e3 * t136 + 0.64e2 * t528 + t532) * t24) * t31) * t112 + t690 * t34 + (0.192e3 * t196 + 0.256e3 * t119 + 0.128e3 * t301 + 0.320e3 * t48) * t157 + (-0.320e3 * t202 + (0.576e3 * t48 + (-0.128e3 * t675 * t8 + t548) * t10) * t14 - 0.576e3 * t5 - 0.64e2 * t8 + (-0.64e2 * t48 * t711 + 0.640e3 * t308) * t10 + (-0.64e2 * t48 * t718 + 0.320e3 * t308) * t15) * t23 + 0.128e3 * t319 + (-t472 + (0.64e2 * t528 * t8 + t272) * t10) * t200 + (t566 + 0.64e2 * t610 * t8 + t327 + (0.128e3 * t681 * t8 - t146 - t659) * t10) * t14 - t608 - 0.384e3 * t149 + 0.128e3 * t743 * t8 + t666 + (0.64e2 * t48 * t753 + 0.384e3 * t257 - t749) * t10 + (0.64e2 * t48 * t764 + 0.64e2 * t257 - t762) * t15 + (t456 * t10 * t317 + (t491 - t490 + (t501 + 0.64e2 * t518 + t532) * t10) * t200 + (-0.704e3 * t136 + 0.64e2 * t568 + t613 + (-0.64e2 * t518 * t7 - t627 + 0.128e3 * t662 + t685) * t10) * t14 + 0.320e3 * t239 - 0.64e2 * t7 + 0.128e3 * t743 + (-0.512e3 + 0.2048e4 * t3 - t252) * t8 + (-0.64e2 * t136 * t792 + 0.64e2 * t5 * t796 - t360 + t801) * t10 + (-0.64e2 * t136 * t806 + 0.64e2 * t810 * t811 - 0.64e2 * t359 + t815) * t15) * t24 + (0.128e3 * t201 + (t10 * t665 - 0.128e3) * t14 - 0.256e3 + t250 + t824 * t10 + t826 * t15) * t158 + (-t388 * t157 + (-0.64e2 * t277 + (-0.128e3 * t615 * t8 + 0.448e3 * t5) * t14 + 0.384e3 * t308 - 0.64e2 * t792 * t48) * t23 + t843 + (0.64e2 * t494 * t8 - t473 - t521) * t200 + (t170 * t643 + 0.128e3 * t619 * t8 + 0.512e3 * t136 - t575) * t14 + 0.640e3 * t257 - t749 + 0.64e2 * t796 * t48 + t824 * t170 + (-t461 * t317 + (t501 - 0.64e2 * t475 + t498) * t200 + (-0.64e2 * t543 * t7 + 0.128e3 * t579 + t624 - t627) * t14 + t360 - 0.64e2 * t711 * t136 + 0.64e2 * t753 * t5 + t801) * t24) * t31 + (-0.64e2 * t387 * t157 + (0.192e3 * t277 + (t272 - t881) * t14 + 0.64e2 * t308 - 0.64e2 * t806 * t48) * t23 - t843 + (-t272 + t881 + t521) * t200 + (t170 * t646 + 0.128e3 * t632 * t8 + t501 - t585) * t14 + 0.320e3 * t257 - t762 + 0.64e2 * t810 * t811 * t8 + t826 * t170 + (t461 * t317 + (-t492 + 0.64e2 * t481 + t507) * t200 + (-0.128e3 * t550 * t7 + 0.128e3 * t589 + t636) * t14 + 0.192e3 * t359 - 0.64e2 * t718 * t136 + 0.64e2 * t764 * t5 + t815) * t24) * t90;
  return(At_fA_re * ((0.256e3 * t12 - 0.256e3 * t17 + (-0.256e3 * t10 * t6 + 0.256e3 * t15 * t6) * t24 + (-0.256e3 * t24 * t6 + 0.256e3 * t26) * t31) * t34 + (-0.384e3 * t12 + 0.384e3 * t17) * t23 + 0.256e3 * t6 * t39 - 0.512e3 * t26 + t51 * t10 - t51 * t15 + (t10 * t64 - t15 * t64 - 0.256e3 * t55 * t6 + 0.512e3 * t6) * t24 + (t24 * t83 + t43 - t72 + t75 - t77) * t31 + (-t24 * t83 - t43 + t72 - t75 + t77) * t90) * PREF_R_CF + Bt_fA_re * t450 * PREF_R_CF + At_fH_re * t919 * PREF_R_CF);
}