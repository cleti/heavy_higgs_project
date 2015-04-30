
#include "AMP_HEADER.h"

double Eval_R_FSR_ISR_SGA (AMP_ARGS, FV const& p3b)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);

 
  double t1;
  double t10;
  double t103;
  double t105;
  double t108;
  double t109;
  double t11;
  double t111;
  double t119;
  double t12;
  double t120;
  double t122;
  double t123;
  double t128;
  double t138;
  double t140;
  double t142;
  double t143;
  double t145;
  double t149;
  double t15;
  double t158;
  double t16;
  double t163;
  double t175;
  double t176;
  double t18;
  double t181;
  double t2;
  double t20;
  double t21;
  double t215;
  double t218;
  double t22;
  double t222;
  double t224;
  double t229;
  double t23;
  double t230;
  double t231;
  double t232;
  double t235;
  double t24;
  double t240;
  double t245;
  double t25;
  double t258;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t34;
  double t35;
  double t36;
  double t37;
  double t38;
  double t39;
  double t4;
  double t40;
  double t41;
  double t42;
  double t43;
  double t44;
  double t45;
  double t46;
  double t48;
  double t49;
  double t5;
  double t51;
  double t53;
  double t55;
  double t56;
  double t57;
  double t59;
  double t6;
  double t63;
  double t64;
  double t65;
  double t67;
  double t68;
  double t69;
  double t7;
  double t71;
  double t73;
  double t74;
  double t75;
  double t76;
  double t78;
  double t8;
  double t80;
  double t83;
  double t84;
  double t85;
  double t87;
  double t9;
  double t90;
  double t92;
  double t93;
  double t94;
  double t96;
  double t98;
  t1 = sp(p1, p2);
  t2 = sp(p1, p3);
  t3 = 0.1e1 / t2;
  t4 = t1 * t3;
  t5 = -t4 + 0.1e1;
  t6 = 0.64e2 * t5;
  t7 = sp(k2, p2);
  t8 = 0.1e1 / t7;
  t9 = t6 * t8;
  t10 = -t6;
  t11 = sp(k2, p1);
  t12 = 0.1e1 / t11;
  t15 = sp(k1, p3);
  t16 = 0.1e1 / t15;
  t18 = sp(k1, p1);
  t20 = sp(p2, p3);
  t21 = 0.1e1 / t20;
  t22 = t1 * t21;
  t23 = t22 - 0.1e1;
  t24 = 0.64e2 * t23;
  t25 = t24 * t8;
  t26 = -t24;
  t30 = sp(k1, p2);
  t32 = sp(p3b, p2);
  t33 = sp(p3b, p3);
  t34 = 0.1e1 / t33;
  t35 = t32 * t34;
  t36 = t35 * t21;
  t37 = sp(p3b, p1);
  t38 = t37 * t34;
  t39 = t38 * t3;
  t40 = -t36 + t39;
  t41 = 0.64e2 * t40;
  t42 = t41 * t1;
  t43 = 0.64e2 * t38;
  t44 = 0.64e2 * t35;
  t45 = t42 - t43 + t44;
  t46 = t45 * t8;
  t48 = -t41 * t1;
  t49 = t48 + t43 - t44;
  t51 = sp(p3b, k1);
  t53 = t51 * t34 * t2;
  t55 = t20 * t51 * t34;
  t56 = -t53 + t55;
  t57 = 0.64e2 * t56;
  t59 = -t57;
  t63 = sp(k2, p3);
  t64 = 0.1e1 / t63;
  t65 = t6 * t64;
  t67 = t24 * t64;
  t68 = t67 * t7;
  t69 = sp(p3b, k2);
  t71 = t69 * t34 * t2;
  t73 = t20 * t69 * t34;
  t74 = -t71 + t73;
  t75 = 0.64e2 * t74;
  t76 = t75 * t64;
  t78 = 0.1e1 / t30;
  t80 = t10 * t64;
  t83 = t26 * t64 * t7;
  t84 = -t75;
  t85 = t84 * t64;
  t87 = 0.1e1 / t18;
  t90 = sp(k1, k2);
  t92 = 0.64e2 * t3;
  t93 = 0.1e1 / t1;
  t94 = 0.64e2 * t93;
  t96 = -0.64e2 * t3 + 0.64e2 * t93;
  t98 = 0.64e2 * t4;
  t103 = t18 * t18;
  t105 = -t96;
  t108 = 0.64e2 * t21;
  t109 = 0.128e3 * t93;
  t111 = 0.64e2 * t21 - 0.64e2 * t93;
  t119 = 0.64e2 * t36;
  t120 = 0.64e2 * t39;
  t122 = -0.64e2 * t38 + 0.64e2 * t35;
  t123 = t122 * t93;
  t128 = t57 * t93;
  t138 = t7 * t11;
  t140 = -t111;
  t142 = t7 * t7;
  t143 = t140 * t64 * t142;
  t145 = -t122 * t93;
  t149 = (t64 * t84 * t93 + t119 - t120 + t145) * t7;
  t158 = t30 * t30;
  t163 = t59 * t93;
  t175 = t11 * t11;
  t176 = t96 * t64 * t175;
  t181 = t75 * t93 * t64;
  t215 = ((t10 * t12 + t9) * t16 * t18 + (t12 * t26 + t25) * t16 * t30 + t46 + t49 * t12 + (t12 * t59 + t57 * t8) * t16 + (t11 * t65 + t42 - t43 + t44 + t68 + t76) * t78 + (t11 * t80 + t43 - t44 + t48 + t83 + t85) * t87) * t90 + (-t92 + t94 + (t7 * t96 + t98 - 0.64e2) * t12) * t16 * t103 + ((t105 * t8 * t11 + t108 + t92 - t109 + t9 + (t111 * t7 - 0.64e2 * t22 + 0.64e2) * t12) * t16 * t30 - t119 + t120 + t123 + ((-t119 + t120 + t123) * t7 + t48 + t43 - t44) * t12 + (t128 + t10 * t8 + (t128 * t7 + 0.64e2 * t53 - 0.64e2 * t55 - t98 + 0.64e2) * t12) * t16 + (t105 * t138 * t64 + t143 + t149) * t78) * t18 + (t11 * t140 * t8 - t108 + t25 + t94) * t16 * t158 + ((t119 - t120 + t145) * t8 * t11 + t119 - t120 + t145 + t46 + (t163 * t8 * t11 + t163 + (-0.64e2 * t22 - 0.64e2 * t53 + 0.64e2 * t55 + 0.64e2) * t8 + t24 * t12) * t16) * t30 + t176 + ((t108 + t92 - t109) * t64 * t7 - t119 + t120 + t123 + t181) * t11 + t143 + t149 + t49 * t8 + t45 * t12 + (t12 * t57 + t59 * t8) * t16 + ((t65 * t7 + t80) * t11 + t67 * t142 + (t42 - t43 + t44 + (-0.64e2 * t22 - 0.64e2 * t71 + 0.64e2 * t73 + 0.64e2) * t64) * t7 + t48 + t43 - t44 + t85) * t78 + ((t176 + (t111 * t64 * t7 - t119 + t120 + t123 + t181) * t11) * t30 + t80 * t175 + (t83 + t48 + t43 - t44 + (-0.64e2 * t4 + 0.64e2 * t71 - 0.64e2 * t73 + 0.64e2) * t64) * t11 + t68 + t42 - t43 + t44 + t76) * t87;
  t218 = 0.128e3 * t5;
  t222 = -t218;
  t224 = 0.128e3 * t23;
  t229 = 0.128e3 * t40;
  t230 = t229 * t1;
  t231 = 0.128e3 * t38;
  t232 = 0.128e3 * t35;
  t235 = 0.128e3 * t56;
  t240 = -t224;
  t245 = -t229 * t1;
  t258 = -0.128e3 * t74;
  return(At_fH_re * t215 * PREF_R_CA + Bt_fA_re * (t218 * t12 * t16 * t103 + ((t12 * t224 + t222 * t8) * t16 * t30 + (t230 - t231 + t232) * t12 + t235 * t12 * t16) * t18 + t240 * t8 * t16 * t158 + ((t245 + t231 - t232) * t8 - t235 * t8 * t16) * t30 + (t222 * t64 * t138 + t240 * t64 * t142 + (t258 * t64 + t231 - t232 + t245) * t7) * t78 + (t218 * t64 * t175 + (t224 * t64 * t7 - t258 * t64 + t230 - t231 + t232) * t11) * t87) * PREF_R_CA);
}
