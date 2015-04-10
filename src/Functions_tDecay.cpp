
#include "../inc/Functions_tDecay.h"






// leading order decay t(S)-> b + (l nu_l)
double Eval_t_blnu(PS_2_3 const& ps, FV& S)
{
  using namespace Constants;
  using namespace RunParameters;
  
  FV const& kl = ps.k1();// lepton
  //FV const& kn = ps.k2();// neutrino
  FV const& kb = ps.k3();// b-quark

  // PRINT_4VEC(kl);
  // PRINT_4VEC(kb);
  // PRINT(Gw4);
  // PRINT(mW2);
  // PRINT(GammaW2);
  
  // set top spin vector in direction of decaying lepton
  // dont forget to adjust kappa if other spin analyzers are used!!!
  S[0] = 0.0;
  S[1] = -kl[1]/kl[0];
  S[2] = -kl[2]/kl[0];
  S[3] = -kl[3]/kl[0];
  // S[0] = 0.0;
  // S[1] = 0.0;
  // S[2] = 0.0;
  // S[3] = 0.0;
  return 4.0*Gw4/(pow(1.0-2.0*kb[0]-mW2,2)+GammaW2*mW2)*(kl[0]*(1.0-2.0*kl[0]));
}

