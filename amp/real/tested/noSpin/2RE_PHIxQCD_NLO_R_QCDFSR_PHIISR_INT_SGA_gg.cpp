double Eval_R_FSR_ISR_SGA (AMP_2_3_ARGS)
{

  AMP_2_3_4VEC_REFS(ps);

  double t1;
  double t10;
  double t100;
  double t101;
  double t102;
  double t103;
  double t104;
  double t105;
  double t106;
  double t109;
  double t11;
  double t110;
  double t111;
  double t112;
  double t113;
  double t116;
  double t118;
  double t119;
  double t12;
  double t120;
  double t122;
  double t124;
  double t127;
  double t129;
  double t13;
  double t131;
  double t132;
  double t133;
  double t135;
  double t141;
  double t15;
  double t150;
  double t153;
  double t16;
  double t162;
  double t163;
  double t17;
  double t171;
  double t173;
  double t174;
  double t177;
  double t181;
  double t182;
  double t184;
  double t185;
  double t186;
  double t187;
  double t19;
  double t2;
  double t20;
  double t207;
  double t21;
  double t210;
  double t211;
  double t214;
  double t216;
  double t219;
  double t22;
  double t23;
  double t25;
  double t26;
  double t27;
  double t272;
  c_double t275;
  double t276;
  double t28;
  double t281;
  double t289;
  double t29;
  double t291;
  double t292;
  double t299;
  double t3;
  double t30;
  double t311;
  double t313;
  double t314;
  double t32;
  double t33;
  double t35;
  double t36;
  double t37;
  double t38;
  double t4;
  double t41;
  double t43;
  double t44;
  double t46;
  double t49;
  double t5;
  double t51;
  double t52;
  double t53;
  double t54;
  double t57;
  double t58;
  double t59;
  double t6;
  double t60;
  double t61;
  double t63;
  double t64;
  double t66;
  double t67;
  double t69;
  double t7;
  double t70;
  double t71;
  double t72;
  double t73;
  double t74;
  double t75;
  double t76;
  double t77;
  double t78;
  double t8;
  double t83;
  double t86;
  double t87;
  double t9;
  double t90;
  double t91;
  double t92;
  double t94;
  double t95;
  double t97;
  t1 = sp(p3, p2);
  t2 = t1 * t1;
  t3 = sp(p3, p1);
  t4 = t1 * t3;
  t5 = -t2 - t4;
  t6 = t1 + t3;
  t7 = 0.1e1 / t6;
  t8 = t5 * t7;
  t9 = 0.1e1 / t3;
  t10 = t8 * t9;
  t11 = 0.1e1 / t1;
  t12 = sp(p3, k1);
  t13 = 0.1e1 / t12;
  t15 = sp(p2, p1);
  t16 = t11 * t13 * t15;
  t17 = t10 * t16;
  t19 = t6 * t6;
  t20 = 0.1e1 / t19;
  t21 = t5 * t20;
  t22 = t21 * t13;
  t23 = 0.128e3 * t22;
  t25 = sp(p2, k2);
  t26 = 0.1e1 / t25;
  t27 = (0.64e2 * t17 - t23) * t26;
  t28 = -t5;
  t29 = t28 * t7;
  t30 = t29 * t9;
  t32 = 0.64e2 * t30 * t16;
  t33 = t28 * t20;
  t35 = 0.128e3 * t33 * t13;
  t36 = t32 - t35;
  t37 = sp(p1, k2);
  t38 = 0.1e1 / t37;
  t41 = sp(p1, k1);
  t43 = sp(p3, k2);
  t44 = 0.1e1 / t43;
  t46 = t11 * t44 * t15;
  t49 = t5 * t44;
  t51 = 0.128e3 * t49 * t20;
  t52 = 0.64e2 * t10 * t46 - t51;
  t53 = sp(p2, k1);
  t54 = 0.1e1 / t53;
  t57 = t3 * t3;
  t58 = t4 + t57;
  t59 = t58 * t7;
  t60 = t59 * t9;
  t61 = t60 * t16;
  t63 = t58 * t20;
  t64 = t63 * t13;
  t66 = 0.64e2 * t61 - 0.128e3 * t64;
  t67 = t66 * t26;
  t69 = t1 - t3;
  t70 = t69 * t7;
  t71 = t9 * t11;
  t72 = t15 * t15;
  t73 = t71 * t72;
  t74 = t70 * t73;
  t75 = 0.64e2 * t74;
  t76 = t69 * t20;
  t77 = t76 * t15;
  t78 = 0.128e3 * t77;
  t83 = t58 * t44;
  t86 = -0.128e3 * t20 * t83 + 0.64e2 * t46 * t60;
  t87 = t86 * t25;
  t90 = -t58;
  t91 = t90 * t7;
  t92 = t91 * t9;
  t94 = 0.64e2 * t92 * t16;
  t95 = t90 * t20;
  t97 = 0.128e3 * t95 * t13;
  t100 = -t69;
  t101 = t100 * t7;
  t102 = t101 * t73;
  t103 = 0.64e2 * t102;
  t104 = t100 * t20;
  t105 = t104 * t15;
  t106 = 0.128e3 * t105;
  t109 = t30 * t46;
  t110 = 0.64e2 * t109;
  t111 = t28 * t44;
  t112 = t111 * t20;
  t113 = 0.128e3 * t112;
  t116 = t92 * t46;
  t118 = t90 * t44;
  t119 = t118 * t20;
  t120 = 0.128e3 * t119;
  t122 = (0.64e2 * t116 - t120) * t25;
  t124 = 0.1e1 / t41;
  t127 = sp(k1, k2);
  t129 = t71 * t13;
  t131 = 0.64e2 * t8 * t129;
  t132 = 0.1e1 / t15;
  t133 = t13 * t132;
  t135 = 0.128e3 * t21 * t133;
  t141 = t41 * t41;
  t150 = t71 * t44;
  t153 = t20 * t132;
  t162 = t2 + 0.2e1 * t4 + t57;
  t163 = t162 * t7;
  t171 = t71 * t15;
  t173 = 0.64e2 * t70 * t171;
  t174 = 0.128e3 * t76;
  t177 = 0.64e2 * t91 * t150;
  t181 = t25 * t25;
  t182 = (-0.128e3 * t118 * t153 + t177) * t181;
  t184 = 0.64e2 * t101 * t171;
  t185 = 0.128e3 * t104;
  t186 = t184 - t185;
  t187 = t186 * t25;
  t207 = 0.64e2 * t8 * t150;
  t210 = -0.128e3 * t153 * t49 + t207;
  t211 = t37 * t37;
  t214 = 0.64e2 * t91 * t129;
  t216 = 0.128e3 * t95 * t133;
  t219 = t53 * t53;
  t272 = ((t36 * t38 + t27) * t41 + t52 * t54 * t37 + t67 * t53 + (t75 - t78) * t26 + (t87 + t75 - t78) * t54 + ((t94 - t97) * t53 + t103 - t106) * t38 + ((t110 - t113) * t37 + t122 + t103 - t106) * t124) * t127 + (t131 - t135 + ((t131 - t135) * t25 + t32 - t35) * t38) * t141 + (((0.64e2 * t129 * t29 - 0.128e3 * t133 * t33) * t26 * t53 + (-0.128e3 * t111 * t153 + 0.64e2 * t150 * t29) * t25 * t54) * t37 + (-0.128e3 * t133 * t162 * t20 + 0.64e2 * t129 * t163 + t27) * t53 + t173 - t174 + t36 * t26 + (t182 + t187) * t54 + (((0.64e2 * t129 * t59 - 0.128e3 * t133 * t63) * t25 + t94 - t97) * t53 + (t173 - t174) * t25 + t103 + (-t185 + t131) * t15 - t23) * t38) * t41 + t210 * t211 + ((t214 - t216) * t26 * t219 + t186 * t26 * t53 + (-0.128e3 * t153 * t162 * t44 + 0.64e2 * t150 * t163) * t25 + t173 - t174 + (t25 * t52 + t110 - t113) * t54) * t37 + (t214 - t216 + t67) * t219 + (t184 - t185 + (t75 + (-t174 + t214) * t15 - t97) * t26) * t53 + t182 + t187 + (t103 - t106) * t26 + (t86 * t181 + (t75 + (-t174 + t177) * t15 - t120) * t25 + t103 - t106) * t54 + (t53 * t66 + t75 - t78) * t38 + ((t210 * t53 + t110 - t113) * t211 + (((0.64e2 * t150 * t59 - 0.128e3 * t153 * t83) * t25 + t173 - t174) * t53 + t122 + t103 + (-t185 + t207) * t15 - t51) * t37 + t87 + t75 - t78) * t124;
  t275 = DenS(k1 + k2, mH, GammaH);
  t276 = RE(t275);
  t281 = 0.128e3 * t17 - 0.256e3 * t22;
  t289 = 0.128e3 * t61 - 0.256e3 * t64;
  t291 = 0.128e3 * t74;
  t292 = 0.256e3 * t77;
  t299 = 0.128e3 * t109 - 0.256e3 * t112;
  t311 = 0.128e3 * t116 - 0.256e3 * t119;
  t313 = 0.128e3 * t102;
  t314 = 0.256e3 * t105;
  return(t272 * FH0 * t276 * PREF_R_CA_At + (t281 * t38 * t141 + (-t281 * t26 * t53 + (t289 * t53 + t291 - t292) * t38) * t41 + t299 * t25 * t54 * t37 - t289 * t26 * t219 + (-t291 + t292) * t26 * t53 + (t311 * t181 + (t313 - t314) * t25) * t54 + (-t299 * t211 + (-t25 * t311 - t313 + t314) * t37) * t124) * FA0 * t276 * PREF_R_CA_Bt);
}
