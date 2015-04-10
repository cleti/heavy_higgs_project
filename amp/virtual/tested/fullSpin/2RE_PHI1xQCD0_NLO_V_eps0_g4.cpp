double Eval_V_2RE_PHI1xQCD0 (
			     PS_2_2 const& ps
			     )
{

#include "EVAL_V_PS_REFS.cpp"

  double t1;
  double t10;
  double t107;
  double t110;
  double t111;
  double t112;
  double t113;
  double t114;
  double t115;
  double t121;
  double t122;
  double t123;
  double t124;
  double t125;
  double t126;
  double t13;
  double t139;
  double t140;
  double t141;
  double t144;
  double t145;
  double t146;
  double t151;
  double t161;
  double t162;
  double t171;
  double t185;
  double t19;
  double t196;
  double t197;
  double t199;
  double t2;
  double t20;
  double t207;
  double t21;
  double t213;
  c_double t22;
  double t223;
  double t225;
  double t23;
  double t24;
  double t248;
  double t250;
  double t253;
  double t254;
  double t255;
  double t261;
  double t267;
  double t28;
  double t282;
  double t287;
  double t29;
  double t3;
  double t30;
  double t303;
  double t34;
  double t35;
  double t36;
  double t37;
  double t4;
  double t40;
  double t41;
  double t42;
  double t43;
  double t47;
  double t5;
  double t53;
  double t55;
  double t56;
  double t57;
  double t6;
  double t61;
  double t62;
  double t68;
  double t69;
  double t7;
  double t71;
  double t77;
  double t8;
  double t81;
  double t83;
  double t9;
  double t96;
  t1 = s * s;
  t2 = CA * t1;
  t3 = At * FH0;
  t4 = 1.0-4.0/s;
  t5 = t3 * t4;
  t6 = Bt * FA0;
  t7 = 0.2e1 * t6;
  t8 = t5 - t7;
  t9 = t2 * t8;
  t10 = beta_y * beta_y;
  t13 = 0.1e1 / (t10 - 0.1e1);
  t19 = VF(0.16e2 * CF * CA * AlphaS3 / 0.3141592653589793e1);
  t20 = t13 * t19;
  t21 = IM(I3_0_0_S12_0_0_0_MU2_0);
  t22 = DenS(s, mH, GammaH);
  t23 = IM(t22);
  t24 = t21 * t23;
  t28 = RE(I3_0_0_S12_0_0_0_MU2_0);
  t29 = RE(t22);
  t30 = t28 * t29;
  t34 = At * FA0;
  t35 = 0.2e1 * t34;
  t36 = Bt * FH0;
  t37 = t35 - t36;
  t40 = RE(I1_MT2_MU2_0);
  t41 = t19 * t40;
  t42 = EPS_(k1, k2, s1, s2);
  t43 = t29 * t42;
  t47 = CF * s;
  t53 = t35 + t36;
  t55 = CF * t53 * t13;
  t56 = sp(s1, k2);
  t57 = t23 * t56;
  t61 = sp(s2, k1);
  t62 = t61 * t19;
  t68 = 0.2e1 * t34 * CA;
  t69 = 0.11e2 / 0.4e1 * t36;
  t71 = (t68 + t69) * t13;
  t77 = t19 * t29;
  t81 = 0.11e2 / 0.4e1 * t5;
  t83 = 0.2e1 * t6 * CA;
  t96 = sp(s2, s1);
  t107 = 0.2e1 * t9 * t20 * t24 - 0.2e1 * t9 * t20 * t30 + 0.4e1 * CF * t37 * t13 * t41 * t43 - 0.2e1 * t47 * t8 * t20 * t40 * t29 + 0.4e1 * t55 * t41 * t57 + 0.4e1 * t55 * t62 * t40 * t23 - 0.8e1 * t71 * t62 * t23 - 0.8e1 * (t68 - t69) * t13 * t77 * t42 + 0.4e1 * s * (t81 - t83) * t20 * t29 - 0.8e1 * t71 * t19 * t23 * t56 - 0.4e1 * s * (t81 + t83) * t13 * t77 * t96 + 0.8e1 * (t83 + 0.11e2 / 0.4e1 * t3) * t13 * t61 * t77 * t56;
  t110 = t3 * t1 * t4 * CF;
  t111 = t4 - 0.1e1;
  t112 = t111 * t13;
  t113 = t112 * t19;
  t114 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t115 = t114 * t23;
  t121 = CF * t111 * s;
  t122 = t3 * t121;
  t123 = t13 * t61;
  t124 = t123 * t19;
  t125 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t126 = t125 * t29;
  t139 = t4 + 0.1e1;
  t140 = CF * t139;
  t141 = t3 + t7;
  t144 = t140 * s * t141 * t13;
  t145 = IM(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t146 = t145 * t23;
  t151 = t34 * t121;
  t161 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t162 = t161 * t29;
  t171 = t19 * t125;
  t185 = 0.2e1 * t110 * t113 * t115 * t96 + 0.4e1 * t122 * t124 * t126 * t56 - 0.2e1 * t110 * t113 * t126 * t96 - 0.4e1 * t122 * t124 * t115 * t56 + 0.2e1 * t144 * t62 * t146 * t56 + 0.8e1 * t151 * t20 * t115 * t42 - 0.2e1 * t110 * t112 * t19 * t114 * t23 - 0.2e1 * t144 * t62 * t162 * t56 - 0.8e1 * t151 * t20 * t126 * t42 + 0.2e1 * t110 * t112 * t171 * t29 - 0.8e1 * t151 * t20 * t125 * t23 * t56 - 0.8e1 * t151 * t123 * t171 * t23;
  t196 = t1 * CF;
  t197 = t5 + t7;
  t199 = t196 * t139 * t197;
  t207 = t140 * s * t37;
  t213 = t140 * s * t53;
  t223 = CA * s;
  t225 = t223 * t141 * t13;
  t248 = -0.8e1 * t151 * t20 * t29 * t114 * t56 - 0.8e1 * t151 * t123 * t77 * t114 + t199 * t20 * t162 * t96 - t199 * t20 * t146 * t96 - 0.2e1 * t207 * t20 * t146 * t42 + 0.2e1 * t213 * t20 * t29 * t145 * t56 + 0.2e1 * t213 * t123 * t77 * t145 + 0.4e1 * t225 * t62 * t24 * t56 + 0.2e1 * t207 * t20 * t162 * t42 + 0.2e1 * t213 * t20 * t161 * t23 * t56 + 0.2e1 * t213 * t123 * t19 * t161 * t23 - 0.4e1 * t225 * t62 * t30 * t56;
  t250 = t196 * t139 * t8;
  t253 = t197 * t13;
  t254 = t2 * t253;
  t255 = t19 * t21;
  t261 = t223 * t37 * t13;
  t267 = t223 * t53 * t13;
  t282 = t19 * t28;
  t287 = t29 * t96;
  t303 = t250 * t20 * t146 - 0.2e1 * t254 * t255 * t23 * t96 - 0.4e1 * t261 * t255 * t23 * t42 + 0.4e1 * t267 * t77 * t21 * t56 + 0.4e1 * t267 * t62 * t29 * t21 - t250 * t20 * t162 + 0.4e1 * t267 * t62 * t28 * t23 + 0.4e1 * t267 * t282 * t57 + 0.2e1 * t47 * t253 * t41 * t287 + 0.4e1 * t261 * t282 * t43 - 0.4e1 * CF * t141 * t123 * t41 * t29 * t56 + 0.2e1 * t254 * t282 * t287;
  return(t107 + t185 + t248 + t303);
} 