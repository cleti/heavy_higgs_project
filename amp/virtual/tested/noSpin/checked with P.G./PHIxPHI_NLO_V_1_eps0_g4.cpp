double Eval_V_PHIxPHI (
		       double const& s, 
		       double const& beta
		       )
{
  double LN_MUR2_S = std::log(MUR2/s);
  double t1;
  double t12;
  double t13;
  double t17;
  double t21;
  double t26;
  double t27;
  double t3;
  double t31;
  double t4;
  double t42;
  double t44;
  double t46;
  double t5;
  double t51;
  double t6;
  double t7;
  double t8;
  t1 = s * s;
  t3 = At * At;
  t4 = t3 * CA;
  t5 = FA0 * FA0;
  t6 = beta * beta;
  t7 = t5 * t6;
  t8 = LN_MUR2_S * LN_MUR2_S;
  t12 = FH0 * FH0;
  t13 = t3 * t12;
  t17 = t5 * Pi2;
  t21 = CA * Pi2;
  t26 = Bt * Bt;
  t27 = t26 * CA;
  t31 = t26 * t12;
  t42 = -CA * t13 * t6 * t8 - CA * t31 * t8 + t13 * t21 * t6 + 0.4e1 * t17 * t4 * t6 - 0.4e1 * t27 * t5 * t8 - 0.4e1 * t4 * t7 * t8 + 0.11e2 * t13 * t6 + 0.4e1 * t17 * t27 + t21 * t31 + 0.16e2 * t27 * t5 + 0.16e2 * t4 * t7 + 0.11e2 * t31;
  t44 = DenS2(s, mH, GammaH);
  t46 = CA * CA;
  t51 = VF(0.8e1 * AlphaS3 * CF * t46 / Pi3);
  return(t1 * s * t42 * t44 * t51 / 0.2e1);

}

