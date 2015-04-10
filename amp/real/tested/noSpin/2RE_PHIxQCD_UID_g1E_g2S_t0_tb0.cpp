double Eval_UID_ES00  (
		      PS_2_3 const& ps,
		      PS_2_2 const& ps_red
		      ) 
{
#include "EVAL_UID_PS_REFS.cpp"

  double t1;
  double t10;
  double t100;
  double t101;
  double t102;
  double t104;
  double t109;
  double t110;
  double t112;
  double t128;
  double t129;
  double t13;
  double t133;
  double t135;
  double t144;
  double t146;
  double t147;
  double t149;
  double t150;
  double t154;
  double t159;
  double t16;
  double t165;
  c_double t175;
  double t176;
  double t179;
  double t19;
  double t2;
  double t20;
  double t21;
  double t22;
  double t23;
  double t25;
  double t26;
  double t3;
  double t33;
  double t35;
  double t36;
  double t37;
  double t38;
  double t4;
  double t42;
  double t44;
  double t45;
  double t46;
  double t47;
  double t5;
  double t50;
  double t54;
  double t56;
  double t57;
  double t6;
  double t63;
  double t68;
  double t69;
  double t7;
  double t70;
  double t72;
  double t74;
  double t8;
  double t82;
  double t9;
  double t93;
  double t94;
  double t96;
  double t98;
  double t99;
  t1 = sp(P1, p3);
  t2 = sp(p2, p1);
  t3 = t1 * t2;
  t4 = sp(K1, p2);
  t5 = t3 * t4;
  t6 = sp(p3, p1);
  t7 = sp(P1, K1);
  t8 = t6 * t7;
  t9 = t8 * t4;
  t10 = sp(K2, p2);
  t13 = sp(K1, p3);
  t16 = sp(K2, p3);
  t19 = sp(p3, p2);
  t20 = t19 * t7;
  t21 = t20 * t2;
  t22 = t4 * t4;
  t23 = t22 * t6;
  t25 = sp(P1, p2);
  t26 = t25 * t6;
  t33 = 0.1e1 / t7;
  t35 = 0.1e1 / t2;
  t36 = 0.1e1 / t4;
  t37 = t35 * t36;
  t38 = EPS_(K1, K2, p2, p3);
  t42 = sp(P1, K2);
  t44 = t19 * t2;
  t45 = t44 * t4;
  t46 = t42 * t1 * t45;
  t47 = t13 * t10;
  t50 = t13 * t25;
  t54 = t26 * t20 * t4;
  t56 = t7 * t2;
  t57 = t56 * t4;
  t63 = t13 * t16;
  t68 = t50 * t21;
  t69 = t16 * t25;
  t70 = t69 * t21;
  t72 = t6 * t19;
  t74 = t42 * t25 * t72 * t4;
  t82 = t1 * t10;
  t93 = 0.2e1 * t1 * t19 * t57;
  t94 = t82 * t21;
  t96 = t1 * t6;
  t98 = t10 * t25 * t96 * t4;
  t99 = t50 * t5;
  t100 = t69 * t5;
  t101 = t46 + 0.2e1 * t47 * t5 - 0.2e1 * t50 * t9 - t54 - 0.2e1 * t13 * t19 * t57 + 0.2e1 * t47 * t26 * t7 - 0.2e1 * t63 * t25 * t7 * t2 - t68 + t70 + t74 - 0.2e1 * t63 * t25 * t2 * t4 + 0.2e1 * t13 * t42 * t45 - 0.2e1 * t82 * t9 - 0.2e1 * t13 * t1 * t57 + 0.2e1 * t16 * t1 * t57 - t93 + t94 + t98 - t99 + t100;
  t102 = t1 * t1;
  t104 = t22 * t102 * t2;
  t109 = t7 * t7;
  t110 = t19 * t19;
  t112 = t109 * t110 * t2;
  t128 = t13 * t13;
  t129 = t128 * t25;
  t133 = t42 * t110 * t56;
  t135 = t22 * t13;
  t144 = t25 * t25;
  t146 = t6 * t4;
  t147 = t13 * t144 * t146;
  t149 = t16 * t144 * t146;
  t150 = t2 * t4;
  t154 = t10 * t102 * t150;
  t159 = t22 * t25 * t96;
  t165 = 0.2e1 * t1 * t22 * t8 - 0.2e1 * t10 * t22 * t96 + 0.2e1 * t16 * t22 * t26 - 0.2e1 * t22 * t42 * t72 + 0.2e1 * t129 * t150 - 0.2e1 * t135 * t26 - 0.2e1 * t135 * t3 + t147 - t149 - t154 - t159;
  t175 = DenS(P1 + p2, mH, GammaH);
  t176 = RE(t175);
  t179 = t104 + t154 - t93 - t94 - t159 - t98 - t99 - t100 - t46 + t112 - t54 - t68 - t70 + t133 + t147 + t149 - t74;
  return((-0.1024e4 * FA0 * (-t10 * t6 * t7 - 0.2e1 * t13 * t2 * t4 - t13 * t2 * t7 + t16 * t2 * t7 - t26 * t4 + t21 + 0.2e1 * t23 + t5 + t9) * t33 * t37 * t38 - 0.256e3 * FH0 * (0.2e1 * t1 * t22 * t4 * t6 + 0.2e1 * t109 * t19 * t4 * t6 - 0.2e1 * t10 * t109 * t72 - 0.2e1 * t109 * t13 * t44 + 0.2e1 * t109 * t16 * t44 + 0.2e1 * t129 * t56 + 0.2e1 * t20 * t23 + t101 + t104 + t112 - t133 + t165) * t33 * t35 * t36) * t176 * PREF_R_CA_At - 0.512e3 * FA0 * t179 * t33 * t37 * t176 * PREF_R_CA_Bt);
}
