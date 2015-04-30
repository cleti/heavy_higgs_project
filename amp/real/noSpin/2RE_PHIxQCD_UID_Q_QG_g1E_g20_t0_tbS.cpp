
#include "UID_HEADER.h"

double Eval_UID_Q_QG_E00S (
			   UID_Q_QG_ARGS)
{

  UID_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);


  double t10;
  double t100;
  double t106;
  double t108;
  double t11;
  double t112;
  double t12;
  double t121;
  double t122;
  double t124;
  double t13;
  double t136;
  double t137;
  double t138;
  double t139;
  double t140;
  double t141;
  double t143;
  double t144;
  double t151;
  double t156;
  double t158;
  double t16;
  double t160;
  double t161;
  double t165;
  double t169;
  double t170;
  double t176;
  double t177;
  double t179;
  double t180;
  double t183;
  double t190;
  double t20;
  double t203;
  double t207;
  double t221;
  double t225;
  double t24;
  double t247;
  double t26;
  double t27;
  double t28;
  double t3;
  double t31;
  double t32;
  double t35;
  double t39;
  double t4;
  double t40;
  double t41;
  double t45;
  double t46;
  double t5;
  double t54;
  double t55;
  double t57;
  double t58;
  double t6;
  double t63;
  double t64;
  double t67;
  double t68;
  double t7;
  double t71;
  double t72;
  double t74;
  double t77;
  double t79;
  double t85;
  double t86;
  double t93;
  double t95;
  t3 = sp(P2, K1);
  t4 = t3 * t3;
  t5 = sp(K2, Q);
  t6 = t4 * t5;
  t7 = sp(K1, Q);
  t10 = Den(P1 + P2, 0);
  t11 = sp(K1, P1);
  t12 = t10 * t11;
  t13 = t12 * Vnd;
  t16 = sp(P2, Q);
  t20 = t4 * t7;
  t24 = t4 * t3;
  t26 = 0.2e1 * t24 * Vdiag;
  t27 = sp(Q, Q);
  t28 = t24 * t27;
  t31 = t5 * t5;
  t32 = t4 * t31;
  t35 = t16 * t16;
  t39 = sp(P2, K2);
  t40 = t39 * t39;
  t41 = t3 * t40;
  t45 = t3 * t7;
  t46 = t5 * t11;
  t54 = t16 * t11;
  t55 = t54 * Vnd;
  t57 = 0.2e1 * t45 * t55;
  t58 = -0.2e1 * Vnd * t39 * t45 * t5 + 0.16e2 * Vdiag * t12 * t41 + 0.2e1 * Vnd * t45 * t46 - 0.24e2 * t13 * t16 * t20 - 0.16e2 * t13 * t16 * t6 + 0.8e1 * t13 * t35 * t4 + 0.8e1 * t13 * t6 * t7 + 0.8e1 * t13 * t28 + 0.8e1 * t13 * t32 - t26 - t57;
  t63 = t3 * t5;
  t64 = t63 * t55;
  t67 = t11 * t39;
  t68 = t67 * Vnd;
  t71 = t7 * t16;
  t72 = sp(K1, K2);
  t74 = t11 * t72 * Vnd;
  t77 = t5 * t16;
  t79 = 0.2e1 * t77 * t68;
  t85 = t3 * t39;
  t86 = t7 * t7;
  t93 = t85 * t5;
  t95 = t11 * Vnd;
  t100 = t16 * t10 * t95;
  t106 = -0.2e1 * t45 * t16 * t39 * Vnd + 0.4e1 * t64 + 0.2e1 * t7 * t5 * t68 + 0.2e1 * t71 * t74 - t79 - 0.2e1 * t77 * t74 - 0.8e1 * t41 * t27 * t13 - 0.8e1 * t85 * t86 * t13 - 0.8e1 * t85 * t35 * t13 - 0.8e1 * t93 * t7 * t10 * t95 + 0.24e2 * t93 * t100 + 0.16e2 * t85 * t7 * t100;
  t108 = t16 * Vnd;
  t112 = t39 * Vnd;
  t121 = t3 * t11;
  t122 = t39 * Vdiag;
  t124 = 0.4e1 * t121 * t122;
  t136 = t27 * t4;
  t137 = t136 * t95;
  t138 = t27 * t3;
  t139 = t40 * Vnd;
  t140 = t138 * t139;
  t141 = -0.16e2 * Vdiag * t10 * t11 * t24 - 0.4e1 * Vdiag * t121 * t72 - 0.2e1 * t11 * t112 * t86 + 0.2e1 * t112 * t3 * t86 - 0.2e1 * t3 * t31 * t95 - 0.2e1 * t108 * t6 + 0.2e1 * t71 * t95 - 0.2e1 * t77 * t95 - t124 - t137 + t140;
  t143 = t27 * t11 * t139;
  t144 = t5 * Vnd;
  t151 = 0.2e1 * t11 * t40 * Vdiag;
  t156 = 0.4e1 * t72 * t4 * Vdiag;
  t158 = 0.4e1 * t85 * Vdiag;
  t160 = 0.2e1 * t41 * Vdiag;
  t161 = t28 * Vnd;
  t165 = t4 * t11 * Vdiag;
  t169 = 0.4e1 * t4 * t39 * Vdiag;
  t170 = 0.4e1 * Vdiag * t67 + 0.2e1 * Vnd * t32 + 0.4e1 * t108 * t20 - 0.2e1 * t144 * t20 + t143 - t151 + t156 - t158 - t160 - t161 + 0.6e1 * t165 + t169;
  t176 = VF(CF * CA2 * AlphaS2);
  t177 = 0.1e1 / t11;
  t179 = 0.1e1 / t3;
  t180 = t176 * t177 * t179;
  t183 = EPS_(K1, K2, P2, Q);
  t190 = t3 * t16;
  t203 = t27 * t72;
  t207 = t72 * t3;
  t221 = t35 * t11;
  t225 = 0.2e1 * Vnd * t203 * t4 + 0.2e1 * Vnd * t207 * t35 - 0.2e1 * Vnd * t221 * t72 - 0.4e1 * Vdiag * t4 - 0.2e1 * t108 * t45 - 0.2e1 * t112 * t136 + 0.2e1 * t138 * t112 - 0.4e1 * t122 * t207 - 0.2e1 * t144 * t190 + t124 - t137 - t140 + t26 + t57 + 0.2e1 * t64 + t79;
  t247 = 0.2e1 * Vnd * t203 * t85 - 0.2e1 * Vnd * t207 * t71 - 0.2e1 * Vnd * t207 * t77 + 0.2e1 * Vnd * t3 * t35 + 0.2e1 * Vnd * t136 - 0.2e1 * Vnd * t221 - 0.2e1 * t138 * t68 + 0.2e1 * t68 * t71 - t143 + t151 - t156 - t158 + t160 - t161 + 0.2e1 * t165 + t169;
  return(-0.8e1 * At_fH_re * (t58 + t106 + t141 + t170) * t180 - 0.32e2 * At_fA_re * t183 * Vnd * (t11 * t7 + 0.4e1 * t12 * t190 + 0.4e1 * t12 * t45 + 0.4e1 * t12 * t63 - t190 - t45 - t46 - t54 + t63) * t176 * t177 * t179 - 0.16e2 * Bt_fA_re * (t225 + t247) * t180);
}
