
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR_INT12 (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);
  AMP_H12_REFS;

  double t1;
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t17;
  double t19;
  double t20;
  double t26;
  double t3;
  double t30;
  double t32;
  double t33;
  double t34;
  double t35;
  double t44;
  double t45;
  double t49;
  double t5;
  double t51;
  double t55;
  double t61;
  double t69;
  c_double t71;
  double t72;
  c_double t74;
  double t75;
  double t77;
  double t79;
  double t9;
  t1 = fA2 * fA1;
  t3 = fH2 * fH1;
  t5 = 0.4096e4 * t1 + 0.1024e4 * t3;
  t9 = t5 * Bt2 * Bt1;
  t11 = sp(p1, p2);
  t12 = t11 * t11;
  t13 = (At1 * At2 * t5 + t9) * t12;
  t14 = sp(k2, p3);
  t15 = 0.1e1 / t14;
  t16 = sp(k1, p3);
  t17 = 0.1e1 / t16;
  t19 = sp(k1, k2);
  t20 = t19 * t19;
  t26 = -0.2048e4 * t1 - 0.512e3 * t3;
  t30 = t26 * Bt2 * Bt1;
  t32 = (At1 * At2 * t26 + t30) * t12;
  t33 = t14 * t14;
  t34 = 0.1e1 / t33;
  t35 = t32 * t34;
  t44 = t16 * t16;
  t45 = 0.1e1 / t44;
  t49 = -t26;
  t51 = t49 * At2 * At1;
  t55 = (Bt1 * Bt2 * t49 + t51) * t12;
  t61 = (t51 + t30) * t12;
  t69 = t13 * t15 * t17 * t20 + (t13 * t15 + t35 + (t13 + (-At1 * At2 * t5 + t9) * t12 * t15) * t17 + t32 * t45) * t19 + (t15 * t55 + t35) * t16 + t13 + t32 * t15 + t61 * t34 + (t14 * t55 + t32) * t17 + (t14 * t32 + t61) * t45;
  t71 = DenS(p1+p2, mH2, GammaH2);
  t72 = RE(t71);
  t74 = DenS(p1+p2, mH1, GammaH1);
  t75 = RE(t74);
  t77 = IM(t71);
  t79 = IM(t74);
  return(PREF_R_PHI_CF * (t69 * t72 * t75 + t69 * t77 * t79));
}
