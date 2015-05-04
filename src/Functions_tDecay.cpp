
#include "../inc/Functions_tDecay.h"


namespace ParametersTDecay
{
  double mW  = 0.0;
  double mW2 = 0.0;
  double GaW  = 0.0;
  double GaW2 = 0.0;  
  double GF  = 0.0;
  double gw2 = 0.0;
  double gw4 = 0.0;
  double GaT = 0.0; // LO top-quark width
  double GaT2= 0.0;
  
  void Init(const HiggsModel& hm)
  {
    using namespace Constants;
    const double& mScale  = hm.Scale();
    const double& mScale2 = hm.Scale2();
    const double& V       = hm.VH();
    const double& mt      = hm.mt();
    const double& mt2     = hm.mt2();
    
    GF  = 1.166E-05*mScale2; // 1.0/(std::sqrt(2.0)*std::pow(V,2));
    mW  = 80.4/mScale;       // 80.385 / mScale;
    mW2 = std::pow(mW,2);
    GaW  = 2.09/mScale;      // 0.208e1 / mScale;
    GaW2 = std::pow(GaW,2);
    gw2  = std::sqrt(32.0)*GF*mW2;
    gw4  = std::pow(gw2,2);
    GaT  = 1.3/mScale;       // GF*std::pow(mt,3)/(std::sqrt(128.0)*Pi)*std::pow(1.0-mW2/mt2,2)*(1.0+2.0*mW2/mt2);
    GaT2 = std::pow(GaT,2);
    std::cout << std::endl << " Init tDecay: " << std::endl;
    PRINT(GF/mScale2);
    PRINT(mW*mScale);
    PRINT(GaW*mScale);
    PRINT(gw2);
    PRINT(GaT*mScale);
  }
}



// leading order decay t(S)-> b + (l nu_l)
double Eval_t_blnu(
		   const PS_2_3& ps,
		   FV& S,
		   const HiggsModel& hm
		   )
{
  using namespace ParametersTDecay;
  static bool INIT = false;
  if (!INIT)
    {
      ParametersTDecay::Init(hm);
      INIT = true;
    }
  
  FV const& kl = ps.k1();// lepton
  //FV const& kn = ps.k2();// neutrino
  FV const& kb = ps.k3();// b-quark

  
  // set top spin vector in direction of decaying lepton
  // dont forget to adjust kappa if other spin analyzers are used!!!
  S[0] = 0.0;
  S[1] = kl[1]/kl[0];
  S[2] = kl[2]/kl[0];
  S[3] = kl[3]/kl[0];
  // S[0] = 0.0;
  // S[1] = 0.0;
  // S[2] = 0.0;
  // S[3] = 0.0;

  // FV Q = kl + kb + kn;
  // PRINT_4VEC(Q);
  // PRINT_4VEC(kb);
  // PRINT_4VEC(kl);
  // PRINT_4VEC(kn);
  // exit(1);
	  
  return gw4/(pow(1.0-2.0*kb[0]-mW2,2)+GaW2*mW2)*(kl[0]*(1.0-2.0*kl[0]))/GaT;
}

