
/*! \file
  \brief Defines the observables used in this program.

  Functions starting with obs_ implement the basic versions of an observables, functions starting with OBS_ implement the interface required for the Distribution/MeanDistribution classes for a specific observable.
  \sa Distribution.h
 */


#ifndef OBSERVABLES_H
#define OBSERVABLES_H


#include "Lorentz.h"

//! Computes invariant mass of two 4-vectors k1 and k2
inline double obs_M12(FV const& k1, FV const& k2)
{
  return sqrt(sp(k1,k1)+sp(k2,k2)+2.0*sp(k1,k2));
}

//! Computes transverse momentum of a 4-vector k1 (transverse = x-y-plane), i.e. k_{T,1} = \sqrt{k1_1}^2+k1_{2}^2}
inline double obs_PT(FV const& k)
{
  return sqrt(k[1]*k[1]+k[2]*k[2]);
}

//! Computes transverse momentum of two 4-vectors k1,k2, k_{T,12} = |k_{T,1}+k_{T,2}|
inline double obs_PT12(FV const& k1, FV const& k2)
{
  return sqrt(pow(k1[1]+k2[1],2)+pow(k1[2]+k2[2],2)); 
}

//! Computes rapidity of a 4-vector k1
inline double obs_Y(FV const& k)
{
  return 0.5*log(((k[0]+k[3])/(k[0]-k[3])));
}

//! Computes rapidity difference of two 4-vectors k1 and k2
inline double obs_DY(FV const& k1, FV const& k2)
{
  return fabs(obs_Y(k1))-fabs(obs_Y(k2));
}

//! Computes the consine of the spatial opening angle of two 4-vectors k1 and k2
inline double obs_PHI(FV const& k1, FV const& k2)
{
  // k1 dot k2 / |k1| / |k2|
  double K1 = LEN(k1);
  double K2 = LEN(k2);
  return (k1[1]*k2[1]+k1[2]*k2[2]+k1[3]*k2[3])/(K1*K2);
}

//! Computes the consine of the opening angle of the spatial projections on the x-y plane
inline double obs_PHIT(FV const& k1, FV const& k2)
{
  // k1_T dot k2_T / |k1_T| / |k2_T|
  double K1 = sqrt(k1[1]*k1[1]+k1[2]*k1[2]);
  double K2 = sqrt(k2[1]*k2[1]+k2[2]*k2[2]);
  return (k1[1]*k2[1]+k1[2]*k2[2])/(K1*K2);
}

//! Computes the spatial triple product of three 4-vectors k1,k2,k3
inline double obs_TriProd(FV const& k1, FV const& k2, FV const& k3)
{
  // = (k1 x k2) dot k3 
  double t3 = k1[1] * k2[2] - k1[2] * k2[1];
  double t7 = -k1[1] * k2[3] + k1[3] * k2[1];
  double t11 = k1[2] * k2[3] - k1[3] * k2[2];
  return (t11 * k3[1] + t3 * k3[3] + t7 * k3[2]);
}

//! Computes the normalized spatial triple product of three 4-vectors k1,k2,k3
inline double obs_TriProdN(FV const& k1, FV const& k2, FV const& k3)
{
  // = (k1 x k2) dot k3 / |k1 x k2| / |k3|
  double t3 = k1[1] * k2[2] - k1[2] * k2[1];
  double t7 = -k1[1] * k2[3] + k1[3] * k2[1];
  double t11 = k1[2] * k2[3] - k1[3] * k2[2];
  double t14 = t3 * t3;
  double t15 = t7 * t7;
  double t16 = t11 * t11;
  double t18 = std::pow(k3[1], 2);
  double t19 = std::pow(k3[2], 2);
  double t20 = std::pow(k3[3], 2);
  double t23 = std::sqrt((t14 + t15 + t16) * (t18 + t19 + t20));
  return (t11 * k3[1] + t3 * k3[3] + t7 * k3[2]) / t23;
}


/*
  Definition of observables used for this project, feel free to add more.
  Two pointers have to be provided as functions arguments:
  a phase space with lab-frame 4-vectors,
  another phase space with parton z.m.f. 4-vectors
*/
class PS_2;

/*!
  \typedef Prototype for all observables. First ps object should contain the lab-frame 4-vectors, the second instance the tt z.m.f. 4-vectors.
*/
typedef double (*OBSFnc)(const PS_2* ps_lab, const PS_2* ps_tt);


//! Computes the top/antitop invariant mass
double OBS_M12   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the top quark transverse momentum (lab frame / parton z.m.f.)
double OBS_PT1   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the antitop transverse momentum (lab frame / parton z.m.f.)
double OBS_PT2   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Copmutes the transverse momentum of top/antitop system (lab frame / parton z.m.f.)
double OBS_PT12  (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the top quark rapidity (lab frame)
double OBS_Y1    (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the antitop quark rapidity (lab frame)
double OBS_Y2    (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the top/antitop rapidity moduli difference Delta Y = |Y_t| - |Y_tbar| (lab frame)
double OBS_DY12  (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the cosine of the opening angle of the top direction and proton beam (tt z.m.f.)
double OBS_Theta1(const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the cosine of the opening angle of the antitop direction and proton beam (tt z.m.f.)
double OBS_Theta2(const PS_2* ps_lab, const PS_2* ps_tt);

/*
  Correlation and polarization observables. They depend on the lepton/antilepton momenta defined in the top/antitop restframes. In this case the matrix elements with full top and antitop spin dependence need to be evaluated.
*/

#ifdef WITH_T_SPIN
//! Computes the opening angle correlation l1.l2 (tt z.m.f.)
double OBS_D12   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the CP odd trilpe correlation |l1 x l2|.k1 (tt z.m.f.)
double OBS_CP1   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the CP odd trilpe correlation |l1 x l2|.k2 (tt z.m.f.)
double OBS_CP2   (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the helicity angle correlation (k1.l1)*(k2.l2) (tt z.m.f.)
double OBS_HEL12 (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the longitudinal polarization k1.l1 (tt z.m.f.)
double OBS_B1    (const PS_2* ps_lab, const PS_2* ps_tt);
//! Computes the longitudinal polarization k2.l2 (tt z.m.f.)
double OBS_B2    (const PS_2* ps_lab, const PS_2* ps_tt);
#endif










#endif
