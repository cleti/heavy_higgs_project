double Eval_V_PHIxPHI (
		       PS_2_2 const& ps
		       )
{

#include "EVAL_V_PS_REFS.cpp"
  
  double t1;
  double t104;
  double t12;
  double t13;
  double t15;
  double t16;
  double t2;
  double t20;
  double t21;
  double t26;
  double t30;
  double t31;
  double t32;
  double t33;
  double t34;
  double t35;
  double t39;
  double t4;
  double t40;
  double t44;
  double t47;
  double t57;
  double t6;
  double t65;
  double t66;
  double t67;
  double t7;
  double t78;
  double t8;
  double t83;
  double t88;
  double t89;
  double t90;
  double t91;
  double t93;
  double t94;
  double t97;
  t1 = At * At;
  t2 = FH0 * FH0;
  t4 = DenS2(s, mH, GammaH);
  t6 = AlphaS3 * CF;
  t7 = CA * CA;
  t8 = 0.1e1 / Pi3;
  t12 = VF(0.8e1 * t6 * t7 * t8);
  t13 = 1.0-4.0/s;
  t15 = s * s;
  t16 = t15 * s;
  t20 = Bt * Bt;
  t21 = t20 * t2;
  t26 = t7 * CA;
  t30 = VF(0.8e1 * t6 * t26 * t8);
  t31 = t4 * t30;
  t32 = t31 * t1;
  t33 = t13 * t2;
  t34 = std::pow(std::log(MUR2/s),2);
  t35 = t16 * t34;
  t39 = FA0 * FA0;
  t40 = t39 * t13;
  t44 = t40 * t16;
  t47 = t31 * t20;
  t57 = t20 * t39 * t16;
  t65 = VF(0.8e1 * t6 * 0.3141592653589793e1 * t26 / Pi2);
  t66 = t4 * t65;
  t67 = t66 * t1;
  t78 = CF * CF;
  t83 = VF(0.8e1 * t78 * t7 * AlphaS3 * t8);
  t88 = t13 * t13;
  t89 = t1 * t88;
  t90 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t91 = t90 * s;
  t93 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t94 = t93 * s;
  t97 = t1 * t13;
  t104 = RE(I1_MT2_MU2_0);
  return(0.11e2 / 0.2e1 * t1 * t2 * t4 * t12 * t13 * t16 + 0.11e2 / 0.2e1 * t21 * t4 * t12 * t16 - t32 * t33 * t35 / 0.2e1 - 0.2e1 * t32 * t40 * t35 + 0.8e1 * t32 * t44 - t47 * t2 * t16 * t34 / 0.2e1 - 0.2e1 * t47 * t39 * t16 * t34 + 0.8e1 * t31 * t57 + 0.2e1 * t67 * t44 + t67 * t33 * t16 / 0.2e1 + 0.2e1 * t66 * t57 + t66 * t21 * t16 / 0.2e1 - t83 * (0.4e1 * t39 + t2) * t16 * t4 * (s * t13 * t20 * t90 + s * t20 * t90 + 0.2e1 * t104 * t20 + 0.2e1 * t104 * t97 + t89 * t91 - 0.2e1 * t89 * t94 + t91 * t97 + 0.2e1 * t94 * t97) / 0.2e1);
}
