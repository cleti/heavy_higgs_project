
#include "AMP_HEADER.h"

double Eval_R_INT_INT (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_R(ap);



  double t1;
  double t10;
  double t103;
  double t109;
  double t110;
  double t128;
  double t13;
  double t130;
  double t131;
  double t132;
  double t133;
  double t139;
  double t14;
  double t148;
  double t153;
  double t156;
  double t16;
  double t17;
  double t173;
  double t175;
  double t176;
  double t177;
  double t185;
  double t186;
  double t188;
  double t189;
  double t194;
  double t2;
  double t206;
  double t212;
  double t213;
  double t215;
  double t216;
  double t22;
  double t228;
  double t23;
  double t238;
  double t239;
  double t247;
  double t25;
  double t258;
  double t26;
  double t269;
  double t271;
  double t273;
  double t282;
  double t283;
  double t295;
  double t3;
  double t306;
  double t32;
  double t33;
  double t359;
  double t36;
  double t365;
  double t37;
  double t38;
  double t382;
  double t386;
  double t39;
  double t394;
  double t4;
  double t40;
  double t41;
  double t42;
  double t43;
  double t44;
  double t46;
  double t47;
  double t5;
  double t58;
  double t6;
  double t64;
  double t65;
  double t67;
  double t73;
  double t78;
  double t8;
  double t82;
  double t83;
  double t84;
  double t85;
  double t86;
  double t87;
  double t88;
  double t9;
  t1 = EPS_(k1, p1, p2, p3);
  t2 = 0.128e3 * t1;
  t3 = sp(p1, p3);
  t4 = t3 * t1;
  t5 = sp(p2, p3);
  t6 = t5 * t1;
  t8 = 0.64e2 * t4 + 0.64e2 * t6;
  t9 = sp(p1, p2);
  t10 = 0.1e1 / t9;
  t13 = sp(k2, p1);
  t14 = 0.1e1 / t13;
  t16 = sp(k1, p2);
  t17 = 0.1e1 / t16;
  t22 = sp(k2, p2);
  t23 = 0.1e1 / t22;
  t25 = sp(k1, p1);
  t26 = 0.1e1 / t25;
  t32 = t9 * t9;
  t33 = 0.1e1 / t32;
  t36 = t3 - t5;
  t37 = t5 + t3;
  t38 = t37 * t37;
  t39 = 0.1e1 / t38;
  t40 = t36 * t39;
  t41 = t40 * t1;
  t42 = 0.2e1 * t5;
  t43 = t42 + t3;
  t44 = 0.1e1 / t37;
  t46 = t1 * t10;
  t47 = t43 * t44 * t46;
  t58 = 0.256e3 * t41;
  t64 = 0.2e1 * t3;
  t65 = t5 + t64;
  t67 = t65 * t44 * t46;
  t73 = t44 * t10;
  t78 = sp(k1, p3);
  t82 = 0.256e3 * t40 * t1 * t9;
  t83 = t3 * t3;
  t84 = t5 * t3;
  t85 = 0.2e1 * t84;
  t86 = t5 * t5;
  t87 = t83 + t85 - t86;
  t88 = t87 * t39;
  t103 = t5 * t44;
  t109 = t83 - t85 - t86;
  t110 = t109 * t39;
  t128 = -0.512e3 * t10 * t103 + 0.256e3 * t33 * t5;
  t130 = t5 * t36;
  t131 = t130 * t39;
  t132 = t86 * t44;
  t133 = t132 * t10;
  t139 = t25 * t25;
  t148 = 0.512e3 * t10 * t3 * t44 - 0.256e3 * t3 * t33;
  t153 = t36 * t44;
  t156 = 0.256e3 * t10 * t153 - 0.512e3 * t40;
  t173 = 0.128e3 * (t83 + 0.6e1 * t84 + t86) * t39 - 0.256e3 * t84 * t73 + 0.256e3 * t84 * t33;
  t175 = t40 * t9;
  t176 = 0.3e1 * t5;
  t177 = -t176 + t3;
  t185 = t9 * t3;
  t186 = t5 * t39;
  t188 = 0.256e3 * t185 * t186;
  t189 = 0.3e1 * t84;
  t194 = -t64 - t42 + t83 + t86;
  t206 = t16 * t16;
  t212 = t3 * t36;
  t213 = t212 * t39;
  t215 = t83 * t44;
  t216 = t215 * t10;
  t228 = t3 + t176;
  t238 = 0.3e1 * t3;
  t239 = t5 + t238;
  t247 = 0.3e1 * t83;
  t258 = 0.3e1 * t86;
  t269 = 0.256e3 * t175;
  t271 = t269 - 0.128e3 * t153;
  t273 = t78 * t78;
  t282 = t32 * t3 * t186;
  t283 = 0.256e3 * t282;
  t295 = t5 * t10;
  t306 = t238 - t5;
  t359 = 0.128e3 * t130 * t10;
  t365 = t39 * t9;
  t382 = 0.128e3 * t212 * t10;
  t386 = 0.256e3 * t40 * t32;
  t394 = 0.512e3 * t282;
  return(Bt_fH_re * ((t10 * t8 - t2) * t14 * t17 + (-t10 * t8 + t2) * t23 * t26) * PREF_R_CA + At_fA_re * ((-0.512e3 * t1 * t14 * t33 + (0.256e3 * t41 - 0.256e3 * t47) * t14 * t17) * t25 + 0.512e3 * t16 * t1 * t23 * t33 + (0.512e3 * t33 * t6 - 0.256e3 * t47 + t58) * t23 + (-0.512e3 * t33 * t4 + t58 + 0.256e3 * t67) * t14 + ((0.512e3 * t4 * t73 - 0.512e3 * t41) * t14 * t78 + (-t82 + 0.256e3 * t88 * t1 - 0.128e3 * (-t64 - t42 + t83 + t85 - t86) * t44 * t46) * t14) * t17 + ((0.256e3 * t41 + 0.256e3 * t67) * t23 * t16 + (-0.512e3 * t103 * t46 - 0.512e3 * t41) * t23 * t78 + (-t82 + 0.256e3 * t110 * t1 - 0.128e3 * (t64 + t42 + t83 - t85 - t86) * t44 * t46) * t23) * t26) * PREF_R_CA + At_fH_re * ((t128 * t14 + (0.128e3 * t131 - 0.128e3 * t133) * t14 * t17) * t139 + ((-t128 * t23 + t14 * t148) * t16 + t156 * t14 * t78 + (-0.256e3 * t33 * t86 + 0.128e3 * t131 + 0.384e3 * t133) * t23 + t173 * t14 + ((-0.128e3 * t177 * t39 * t5 + 0.128e3 * t130 * t73 - 0.128e3 * t175) * t14 * t78 + (-t188 + 0.128e3 * t5 * (t64 + t42 + t83 + t189) * t39 + 0.64e2 * t5 * t194 * t73) * t14) * t17) * t25 - t148 * t23 * t206 + (-t156 * t23 * t78 + t173 * t23 + (-0.256e3 * t33 * t83 - 0.128e3 * t213 + 0.384e3 * t216) * t14) * t16 + ((-0.128e3 * t175 - 0.128e3 * t5 * (t3 + 0.5e1 * t5) * t39 + 0.128e3 * t5 * t228 * t73) * t23 + (0.128e3 * t175 - 0.128e3 * t3 * (0.5e1 * t3 + t5) * t39 + 0.128e3 * t3 * t239 * t73) * t14) * t78 + (-t188 + 0.128e3 * t5 * (-t64 - t42 + t247 + t84) * t39 - 0.64e2 * t5 * (-t64 - t42 + t247 - t86) * t73) * t23 + (-t188 + 0.128e3 * t3 * (-t64 - t42 + t84 + t258) * t39 + 0.64e2 * t3 * (t64 + t42 + t83 - t258) * t73) * t14 + (t271 * t14 * t273 + (t269 - 0.64e2 * (t64 - t42 + t83 + t85 + t86) * t44) * t14 * t78 + (t283 - 0.256e3 * t84 * (0.2e1 + t3 + t42) * t39 * t9 + 0.128e3 * t84 * (0.2e1 + t176 + t3) * t44 - 0.128e3 * (0.2e1 + t5) * t3 * t295) * t14) * t17 + ((-0.128e3 * t213 - 0.128e3 * t216) * t23 * t206 + ((0.128e3 * t3 * t306 * t39 - 0.128e3 * t212 * t73 + 0.128e3 * t175) * t23 * t78 + (-t188 + 0.128e3 * t3 * (t64 + t42 + t189 + t86) * t39 + 0.64e2 * t3 * t194 * t73) * t23) * t16 - t271 * t23 * t273 + (-t269 - 0.64e2 * (-t64 + t42 + t83 + t85 + t86) * t44) * t23 * t78 + (t283 - 0.256e3 * t84 * (0.2e1 + t64 + t5) * t39 * t9 + 0.128e3 * t84 * (0.2e1 + t5 + t238) * t44 - 0.128e3 * (0.2e1 + t3) * t3 * t295) * t23) * t26) * PREF_R_CA + Bt_fA_re * ((0.256e3 * t44 * t5 * t9 - 0.256e3 * t132 - t359) * t14 * t17 * t25 + (0.256e3 * t5 * t306 * t365 - 0.256e3 * t5 * (t64 - t5) * t44 + t359) * t23 + (-0.256e3 * t3 * t177 * t365 + 0.256e3 * t3 * (t3 - t42) * t44 - t382) * t14 + ((0.128e3 * t44 * t87 - 0.256e3 * t88 * t9 + t386) * t14 * t78 + (0.512e3 * t39 * t43 * t84 * t9 + 0.256e3 * t10 * t3 * t86 - 0.256e3 * t228 * t44 * t84 - t394) * t14) * t17 + ((0.256e3 * t185 * t44 - 0.256e3 * t215 + t382) * t23 * t16 + (-0.128e3 * t109 * t44 + 0.256e3 * t110 * t9 - t386) * t23 * t78 + (0.512e3 * t39 * t65 * t84 * t9 + 0.256e3 * t10 * t5 * t83 - 0.256e3 * t239 * t44 * t84 - t394) * t23) * t26) * PREF_R_CA);
}
