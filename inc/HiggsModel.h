
/*! \file
 \brief  Model specific settings
*/ 

#ifndef HIGGSMODEL_H
#define HIGGSMODEL_H

#include <iterator>
#include <memory>
#include <complex>
#include <iostream>
#include <iomanip>

#include "ScalarIntegrals.h"




/* class num_except : std::exception { */
/* public: */
/*   std::string location = ""; */
/*   std::string varname = ""; */
/*   double val = 0.0; */
/*   const char* what() const noexcept */
/*   { */
/*     char* buffer = new char[100]; */
/*     snprintf(buffer,100,"Numeric exception in %s. Value of %s = %e not permitted.",location.c_str(),varname.c_str(),val); */
/*     return buffer; */
/*   } */
/* }; */



/*!
\typedef alias for complex values
 */
using c_double = std::complex<double>;

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/*!
  This class contains the parameters that describe a single neutral Higgs boson.
  \sa HiggsModel
*/
class HiggsBoson
{
 private:
  //! mass
  double d_M;
  //! mass squared
  double d_M2;
  //! width
  double d_G;
  //! width squared
  double d_G2;

  //! reduced scalar-top coupling diveded by combined vacuum expectation value
  double d_At;
  //! reduced pseudoscalar-top coupling diveded by combined vacuum expectation value
  double d_Bt;
  //! reduced scalar-bottom coupling diveded by combined vacuum expectation value
  double d_Ab;
  //! reduced pseudoscalar-bottom coupling diveded by combined vacuum expectation value
  double d_Bb;
  
  //! effective coupling to gluons
  c_double d_FH_eff;
  c_double d_FA_eff;

  //! current value of c.m.e. used to compute the form factors
  double   d_S_FF;
  //! gg-scalar 1-loop form factors
  /*!
    \sa SetFormFactors
  */
  c_double d_FH_full;
  //! gg-pesudoscalar 1-loop form factors
  /*!
    \sa SetFormFactors
  */
  c_double d_FA_full;

  
  //! current value of c.m.e. used to compute the propagator
  double   d_S_Den;
  //! Higgs propagator squared
  /*!
    \sa SetPropagator
  */
  double   d_DenSq;
  //! Higgs propagator
  /*!
    \sa SetPropagator
  */
  c_double d_Den;
  
public:
  HiggsBoson(double const& M,
	     double const& G,
	     double const& Vh,
	     double const& At,
	     double const& Bt,
	     double const& Ab=0.0,
	     double const& Bb=0.0);
  ~HiggsBoson() {}
  
  double const& M()  const { return d_M; }
  double const& G()  const { return d_G; }
  double const& At() const { return d_At; }
  double const& Bt() const { return d_Bt; }  
  double const& Ab() const { return d_Ab; }
  double const& Bb() const { return d_Bb; }
  
  c_double const& GetFH0() const { return d_FH_eff; }
  c_double const& GetFA0() const { return d_FA_eff; }

  //! Compute 1-loop form factors for given c.m.e., top- and bottom mass and store values in member variables. Note that the form factors are recomputed only if the c.m.e. changes (compared to last call)! The values stored in the respective member variables are - F_s / (4 s) and F_p / (8 s).
  /*!
    \param S c.m.e.
    \param mt2 top-mass squared
    \param mb2 bottom-mass squared
  */
  void SetFormFactors(double const& S,double const& mt2,double const& mb2);
  c_double const& GetFs() const { return d_FH_full; }
  c_double const& GetFp() const { return d_FA_full; }

  //! returns effective gg-scalar coupling if EFF=true, full 1-loop form factor otherwise
  c_double const& GetFH(bool EFF) const { if (EFF) { return d_FH_eff; } else { return d_FH_full; } }
  //! returns effective gg-pseudoscalar coupling if EFF=true, full 1-loop form factor otherwise
  c_double const& GetFA(bool EFF) const { if (EFF) { return d_FA_eff; } else { return d_FA_full; } }

  //! compute propagator value for given c.m.e. and store values in member variables
    /*!
    \param S momentum squared in the Higgs propagators
  */
  void SetPropagator(double const& S);
  double   const& GetPropagatorSq() const { return d_DenSq; }
  c_double const& GetPropagator()   const { return d_Den; }


};
/*!
\typedef HiggsBoson pointer
 */
using HPtr = std::shared_ptr<HiggsBoson>;
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/*!
This structure contains the Higgs specific prefactors, i.e. couplings and the propagator denominator. The prefactors depend on the momentum of the Higgs boson and have to be reset whenever the phase space point changes (usually for every call of the integrand). This is done via the HiggsModel class.
\sa HiggsModel::SetPrefactors()
*/
struct HiggsPrefactors {
 public:
  // PHIxQCD amplitudes //////
  double At_fH_re = 0.0;		
  double At_fA_re = 0.0;		
  double Bt_fH_re = 0.0;
  double Bt_fA_re = 0.0;
  double At_fH_im = 0.0;
  double At_fA_im = 0.0;
  double Bt_fH_im = 0.0;
  double Bt_fA_im = 0.0;
  // PHIxPHI amplitudes //////
  double At2_fH2_De = 0.0;
  double At2_fA2_De = 0.0;
  double Bt2_fH2_De = 0.0;
  double Bt2_fA2_De = 0.0;
  double At_Bt_fH2_De = 0.0;
  double At_Bt_fA2_De = 0.0;
  double At_Bt_fH2_DeIM = 0.0;
  double At_Bt_fA2_DeIM = 0.0;

  void Reset();
  void Print(std::ostream& ost);
};
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/*!
This structure contains the prefactors used in the amplitudes. The prefactors depend on AlphaS and have to be reset whenever it changes (usually only once for each run). This is done via the HiggsModel class.
\sa HiggsModel::SetAmpPrefactors()
*/
struct AmplitudePrefactors {
 public:
  /* born */
  double PREF_B_PHIxPHI = 0.0;
  double PREF_B_PHIxQCD = 0.0;
  double PREF_B_QCDxQCD = 0.0;
  double PREF_B_QCDxQCD_CF = 0.0;
  double PREF_B_QCDxQCD_CA = 0.0;
  double PREF_B_QCDxQCD_CFCA = 0.0;
  /* virtual */
  double PREF_V_PHIxQCD    = 0.0;
  double PREF_V_PHIxQCD_CF = 0.0;
  double PREF_V_PHIxQCD_CA = 0.0;
  double PREF_V_PHIxQCD_CFCA2 = 0.0;
  double PREF_V_PHIxQCD_Nf = 0.0;
  double PREF_V_PHIxQCD_CT = 0.0;
  double PREF_V_PHIxPHI = 0.0;
  double PREF_V_PHIxPHI_CA = 0.0;
  double PREF_V_PHIxPHI_CF = 0.0;
  /* real */
  double PREF_R_PHIxQCD = 0.0;
  double PREF_R_PHIxQCD_CF = 0.0;
  double PREF_R_PHIxQCD_CA = 0.0;
  double PREF_R_PHIxQCD_CFCA2 = 0.0;
  double PREF_R_PHIxPHI = 0.0;
  double PREF_R_PHIxPHI_CA = 0.0;
  double PREF_R_PHIxPHI_CF = 0.0;
  /* UID */
  double PREF_UID_TF = 0.0;
  double PREF_UID_CA = 0.0;
  double PREF_UID_CF = 0.0;

  /* Needed in  EVAL_V  */
  double MUR2   = 0.0;
  double AlphaS = 0.0;
  
  void Print(std::ostream& ost);
};
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/*!
This class contains all the physical, model specific parameters, i.e. the strong coupling AlphaS, renormalization and factorization scales MUR, MUF, the third generation quark masses mt and mb as well as the combined Higgs vacuum expectation value VH and the individual Higgs boson parameters. It also provides appropriate setter functions. It contains instances of the HiggsPrefactors and Amprefactors structures that are needed for the evaluation of the amplitudes. Take care to provide numerical values consistently in the same units.
\sa AmpPrefactors, HiggsPrefactors, HiggsBoson
*/
class HiggsModel
{
 private:
  std::string d_Name;
  //! strong couplin
  double d_AlphaS;
  //! strong coupling squared
  double d_AlphaS2; 
  
  //! renormalization scale
  double d_MUR;  
  //! renormalization scale square
  double d_MUR2;
  //! factorization scale
  double d_MUF; 
  //! factorization scale squared
  double d_MUF2; 

  //! top-quark mass
  double d_mt; 
  //! top-quark mass squared
  double d_mt2; 
  //! bottom-quark mass
  double d_mb;  
  //! bottom-quark mass squared
  double d_mb2; 
  //! combined Higgs VEV
  double d_VH;
  //! mass scale for normalization of dimensionful quantities
  double d_Scale;
  double d_Scale2;
  
  //! contains the pointers to the Higgs bosons
  std::vector<HPtr> d_Bosons;
  //! Higgs specific prefactors (couplings, propagator)
  HiggsPrefactors     d_HPref; 
  //! Amplitude prefactors (AlphaS,Pi,CF,CA,...)
  AmplitudePrefactors d_APref; 

 public:
  explicit HiggsModel(std::string const & name="noname");
  ~HiggsModel();

  std::string const& Name() const { return d_Name; }
  double const& AlphaS()    const { return d_AlphaS; }
  double const& AlphaS2()   const { return d_AlphaS2; }
  double const& MUR() 	    const { return d_MUR; }
  double const& MUR2()	    const { return d_MUR2; }
  double const& MUF() 	    const { return d_MUF; }
  double const& MUF2()	    const { return d_MUF2; }
  double const& mt()  	    const { return d_mt; }
  double const& mt2() 	    const { return d_mt2; }
  double const& mb()  	    const { return d_mb; }
  double const& mb2() 	    const { return d_mb2; }
  double const& VH()  	    const { return d_VH; }
  double const& Scale()     const { return d_Scale; }
  double const& Scale2()    const { return d_Scale2; }  
  int NBosons()             const { return d_Bosons.size(); }
  
  //! Use this member to change AlphaS. It automatically resets the values of the amplitude prefactors.
  void SetAlphaS(double const& val) { d_AlphaS = val; d_AlphaS2 = std::pow(val,2); SetAmpPrefactors(); }
  void SetMUR(double const& val)    { d_MUR = val/d_Scale; d_MUR2 = std::pow(d_MUR,2); d_APref.MUR2 = d_MUR2; }
  void SetMUF(double const& val)    { d_MUF = val/d_Scale; d_MUF2 = std::pow(d_MUF,2); }
  void SetMt(double const& val)     { d_mt  = val/d_Scale; d_mt2  = std::pow(d_mt,2); }
  void SetMb(double const& val)     { d_mb  = val/d_Scale; d_mb2  = std::pow(d_mb,2); }
  void SetVH(double const& val)     { d_VH  = val/d_Scale; }
  void SetScale(double const& val)  { d_Scale  = val; d_Scale2 = val*val; }
  
  //! Reset values of the Higgs prefactors. All Higgs bosons in the vector d_Bosons will be taken into account.
  /*!
    \param S momentum squared in the Higgs propagators
    \param EFF use effective Higgs-top coupling (large mt limit) if true. Couplings to the bottom-quark have no effect in this case. The full one-loop form factors are used otherwise.
  */
  void SetHiggsPrefactors(double const& S, bool EFF);
  HiggsPrefactors     const& GetHiggsPrefactors() const { return d_HPref; }

  //! Reset values of the amplitude prefactors. 
  void SetAmpPrefactors();
  AmplitudePrefactors const& GetAmpPrefactors()   const { return d_APref; }

  //! Add a Higgs boson to the vector d_Bosons. The dimensionful parameters M and G will be rescaled with d_Scale.
  /*!
    \param M mass
    \param G width
    \param a_t reduced scalar coupling to the top-quark
    \param b_t reduced pseudoscalar coupling to the top-quark
    \param a_t reduced scalar coupling to the bottom-quark
    \param b_t reduced pseudoscalar coupling to the bottom-quark
  */
  void AddBoson(
		double const& M,
		double const& G,
		double const& a_t=1.0,
		double const& b_t=1.0,
		double const& a_b=0.0,
		double const& b_b=0.0);

  //! Add a scalar Higgs boson to the vector d_Bosons. The dimensionful parameters M and G will be rescaled with d_Scale.
  /*!
    \param M mass
    \param G width
    \param a_t reduced scalar coupling to the top-quark
    \param a_t reduced scalar coupling to the bottom-quark
  */  
  void AddScalar(
		 double const& M,
		 double const& G,
		 double const& a_t=1.0,
		 double const& a_b=0.0);

  //! Add a pseudoscalar Higgs boson to the vector d_Bosons. The dimensionful parameters M and G will be rescaled with d_Scale.
  /*!
    \param M mass
    \param G width
    \param b_t reduced pseudoscalar coupling to the top-quark
    \param b_t reduced pseudoscalar coupling to the bottom-quark
  */    
  void AddPseudoscalar(
		       double const& M,
		       double const& G,
		       double const& b_t=1.0,
		       double const& b_b=0.0);

  //! Remove the last Higgs boson added to the vector d_Bosons.
  void PopBoson();
  
  //! Print the model parameters.
  /*!
    \param ost output stream
    \param mScale mass scale, used to restore proper normalization of quantities with mass dimension
  */   
  void Print(std::ostream& ost, double const& mScale);
};
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////



#endif
