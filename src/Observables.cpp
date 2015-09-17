


#include "../inc/Observables.h"
#include "../inc/PhaseSpace.h"




/////////////////////////////////////////////////////////////////////////////////////////////
// definition of observables used for this project, feel free to add more ///////////////////
// no inlining here, since the functions need to get addressed            ///////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
double OBS_M12(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_M12(ps_lab->k1(),ps_lab->k2());
}

double OBS_PT1(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_PT(ps_lab->k1());
}

double OBS_PT2(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_PT(ps_lab->k2());
}

double OBS_PT12(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_PT12(ps_lab->k1(),ps_lab->k2());
}

double OBS_Y1(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_Y(ps_lab->k1());
}
double OBS_Y2(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_lab == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  return obs_Y(ps_lab->k2());
}

double OBS_DY12(const PS_2* ps_lab, const PS_2* ps_tt)
{
  return obs_DY(ps_lab->k1(),ps_lab->k2());
}

// cosine of Collin-Soper angle defined by top quark and proton beam direction in tt z.m.f.
double OBS_Theta1(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // in case of 2->n reactions with n>2, the tt system can have nonzero transverse momentum
  // we have to compute the vector that dissects the angle between P1 and -P2 (proton momenta)
  if (ps_tt->whattype()>22)
    {
      // the sum of the normalized vectors dissects the angle between them
      FV N = ps_tt->P1().norm3() - ps_tt->P2().norm3();
      
      if (LEN(N)<1e-10)
      	{
      	  WARNING("LEN(N)<1e-10");
      	}
      double cphi1 = obs_PHI(N, ps_tt->P1());
      double cphi2 = obs_PHI(N,-ps_tt->P2());
      if (std::fabs(cphi1-cphi2)>1e-12)
	{
	  PRINT_4VEC(N);
	  PRINT(cphi1);
	  PRINT(cphi2);
	  PRINT(obs_PHI(N,ps_tt->k1()));
	  PRINT(obs_PHI(N,ps_tt->k2()));
	  FV K = ps_tt->k1()+ps_tt->k2();
	  FV Q = ps_tt->P1()+ps_tt->P2();
	  PRINT_4VEC(ps_tt->k1());
	  PRINT_4VEC(ps_tt->k2());
	  PRINT_4VEC(ps_tt->P1().norm3());
	  PRINT_4VEC(ps_tt->P2().norm3());      
	  PRINT_4VEC(K);
	  PRINT_4VEC(Q);
	  EXIT(1);
	}
      return obs_PHI(ps_tt->k1(),N);
    }
  return obs_PHI(ps_tt->k1(),ps_tt->P1());
}

// cosine of Collin-Soper angle defined by antitoptop quark and proton beam direction in tt z.m.f.
double OBS_Theta2(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // in case of 2->3, where the tt system can have transverse momentum
  // we have to compute the vector that dissects P1 and -P2 (proton momenta)
  if (ps_tt->whattype()>22)
    {
      FV N = ps_tt->P1().norm3() - ps_tt->P2().norm3();
      // PRINT_4VEC(N);
      // PRINT(obs_PHI(N,ps_tt->P1()));
      // PRINT(obs_PHI(N,ps_tt->P2()));
      // FV K = ps_tt->k1()+ps_tt->k2();
      // FV Q = ps_tt->P1()+ps_tt->P2();
      // PRINT_4VEC(K);
      // PRINT_4VEC(Q);
      // EXIT(1);
      return obs_PHI(ps_tt->k2(),N);
    }
  return obs_PHI(ps_tt->k2(),ps_tt->P1());
}




/*
  Computation with full spin dependence:
  in this case it is assumed that the 4-vectors s1_r and s2_r in each ps instance
  contain the lepton/antilepton momenta in the top/antitop restframe!
*/
#ifdef WITH_T_SPIN  
double OBS_D12(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  return (-3.0)*obs_PHI(ps_tt->s1_r(),ps_tt->s2_r())/Constants::kappa_p/Constants::kappa_m;
}

double OBS_CP1(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1 is the top momentum in the ttbar zero-momentum frame
  return obs_TriProdN(ps_tt->s1_r(),ps_tt->s2_r(),ps_tt->k1())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_CP2(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1 is the top momentum in the ttbar zero-momentum frame
  return obs_TriProdN(ps_tt->s1_r(),ps_tt->s2_r(),ps_tt->k2())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_HEL12(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (-9.0)*obs_PHI(ps_tt->k1(),ps_tt->s1_r())*obs_PHI(ps_tt->k2(),ps_tt->s2_r())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_B1(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (3.0)*obs_PHI(ps_tt->k1(),ps_tt->s1_r())/Constants::kappa_p;
}
double OBS_B2(const PS_2* ps_lab, const PS_2* ps_tt)
{
#ifdef DEBUG
  if (ps_tt == nullptr)
    {
      ERROR_NOEXIT("nullptr received in observable");
      return 0.0;
    }
#endif
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (3.0)*obs_PHI(ps_tt->k2(),ps_tt->s2_r())/Constants::kappa_m;
}
#endif
