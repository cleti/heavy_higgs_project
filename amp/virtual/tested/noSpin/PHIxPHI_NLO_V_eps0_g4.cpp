
#include "../../../../inc/Functions_Shared.h"

double Eval_V_PHIxPHI (
		       PS_2_2 const& ps
		       )
{

#include "EVAL_V_PS_REFS"
  
  double t1;
  double t12;
  double t13;
  double t17;
  double t2;
  double t21;
  double t26;
  double t27;
  double t3;
  double t31;
  double t4;
  double t42;
  double t44;
  double t46;
  double t47;
  double t5;
  double t51;
  double t54;
  double t59;
  double t6;
  double t64;
  double t65;
  double t66;
  double t67;
  double t69;
  double t7;
  double t70;
  double t73;
  double t8;
  double t80;
  t1 = s * s;
  t2 = t1 * s;
  t3 = At * At;
  t4 = t3 * CA;
  t5 = FA0 * FA0;
  t6 = 1.0-4.0/s;
  t7 = t5 * t6;
  t8 = std::pow(std::log(MUR2/s),2);
  t12 = FH0 * FH0;
  t13 = t3 * t12;
  t17 = t5 * Pi2;
  t21 = CA * Pi2;
  t26 = Bt * Bt;
  t27 = CA * t26;
  t31 = t26 * t12;
  t42 = -CA * t13 * t6 * t8 - CA * t31 * t8 + t13 * t21 * t6 + 0.4e1 * t17 * t4 * t6 - 0.4e1 * t27 * t5 * t8 - 0.4e1 * t4 * t7 * t8 + 0.11e2 * t13 * t6 + 0.4e1 * t17 * t27 + t21 * t31 + 0.16e2 * t27 * t5 + 0.16e2 * t4 * t7 + 0.11e2 * t31;
  t44 = DenS2(s, mH, GammaH);
  t46 = CA * CA;
  t47 = 0.1e1 / Pi3;
  t51 = VF(0.8e1 * AlphaS3 * CF * t46 * t47);
  t54 = CF * CF;
  t59 = VF(0.8e1 * t54 * t46 * AlphaS3 * t47);
  t64 = t6 * t6;
  t65 = t3 * t64;
  t66 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0);
  t67 = t66 * s;
  t69 = RE(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1);
  t70 = t69 * s;
  t73 = t3 * t6;
  t80 = RE(I1_MT2_MU2_0);
  return(t2 * t42 * t44 * t51 / 0.2e1 - t59 * (0.4e1 * t5 + t12) * t2 * t44 * (s * t26 * t6 * t66 + s * t26 * t66 + 0.2e1 * t26 * t80 + t65 * t67 - 0.2e1 * t65 * t70 + t67 * t73 + 0.2e1 * t70 * t73 + 0.2e1 * t73 * t80) / 0.2e1);
}

