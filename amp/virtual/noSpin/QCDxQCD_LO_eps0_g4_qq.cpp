
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_QQ (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_B(ap);


  double t1;
  double t2;
  t1 = beta * beta;
  t2 = y * y;
  // -----|-- sum over top spins
  return 4.0*((0.8e1 * t1 * t2 - 0.8e1 * t1 + 0.16e2) * PREF_B_QCDxQCD);
}
