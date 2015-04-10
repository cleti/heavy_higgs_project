double Eval_B_QCDxQCD_QQ (
		       PS_2_2 const& ps
		       )
{

#include "EVAL_V_PS_REFS"

  double t1;
  double t2;
  t1 = 0.1e1 - 0.4e1 / s;
  t2 = beta_y * beta_y;
  return(0.8e1 * PREF_B_QCDxQCD * (0.2e1 - (t1 - t2)));
}


