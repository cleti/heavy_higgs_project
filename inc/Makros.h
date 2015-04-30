


#ifndef MAKROS_H
#define MAKROS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <boost/utility.hpp>

#define CHECKNAN(var)  if (var != var) {std::cout << std::endl << " !!! in " << __FUNCTION__  << ", line " <<  __LINE__ << std::endl;ps.print();std::cout << std::endl << " !!! " << #var << " is nan !!! " << std::endl; std::exit(1);}
#define CHECKNA(var) if ( std::isnan(var) || std::isinf(var) ) {std::cout << std::endl << " !!! in " << __FUNCTION__  << ", line " <<  __LINE__ << std::endl << std::endl << " !!! " << #var << " is nan !!! " << std::endl; std::exit(1);}


#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)



#define WARNING(MSG)     std::cout << std::endl << " Warning in " << __FUNCTION__  << ", line " <<  __LINE__ << ": " <<  MSG << ". " << std::endl;


#define ERROR(MSG)       std::cout << std::endl << " ! ERROR in " << __FUNCTION__  << ", line " <<  __LINE__ << ": " <<  MSG << "! " << std::endl; std::exit(1);

#define SLEEP(S)         std::this_thread::sleep_for( std::chrono::seconds(S) ); 

#define EXIT(S) std::cout << std::endl << " exit called from " << __FUNCTION__ << " line " <<   __LINE__ << std::endl; exit(S);

#define ERRORN(NAME,MSG) std::cout << std::endl << " ! ERROR in " << NAME << ", " << __FUNCTION__  << ", line " <<  __LINE__ << ": " <<  MSG << "! " << std::endl; std::exit(1);
#define PRINT(VAR)       std::cout << std::endl << " Value of " << std::setw(20) <<  #VAR << " = " << std::setw(20)   << VAR << std::endl;
#define PRINTS(OST,VAR)       OST << std::endl << " Value of " << std::setw(20) <<  #VAR << " = " << std::setw(20)   << VAR << std::endl;

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define LOG(MSG) MyLog::gLogger().log(MSG "\n",0);
#define LOG_AT(MSG) MyLog::gLogger().log(AT,0);LOG(MSG)
#define LOG_ONCE(MSG) { static int msg = 1; if (unlikely(msg)) { LOG_AT(MSG); msg = 0; } }
#define LOG_VAL(VAR) {std::stringstream O; O << "Value of " << #VAR << " = " << VAR; MyLog::gLogger().log(O.str());}

#define WARN_ONCE(MSG) { static int msg = 1; if (unlikely(msg)) { WARNING(MSG); msg = 0; } }

#define WARN(fnc,msg) printf("\n Warning in %s : %s line %i, %s. \n",__FILE__,#fnc,__LINE__,msg); 

#define RE(X) X.real()
#define IM(X) X.imag()
#define VF(X) (X)

// Lorentz product of two 4 vectors K1,K2
#define sp(K1,K2) (K1[0]*K2[0]-(K1[1]*K2[1]+K1[2]*K2[2]+K1[3]*K2[3]))
// return mass squared of a 4-vector
#define MSQ(K) sp(K,K)
#define LEN(K) sqrt(K[1]*K[1]+K[2]*K[2]+K[3]*K[3])


#define PRINT_4VEC(P)							\
  {									\
    std::cout << std::endl << #P << " =  [ ";				\
      std::cout << std::setw(8) << P[0] << ", ";			\
      std::cout << std::setw(8) << P[1] << ", ";			\
      std::cout << std::setw(8) << P[2] << ", ";			\
      std::cout << std::setw(8) << P[3] << " ]";			\
      std::cout << std::endl << #P << "^2 = " << std::setw(8) << MSQ(P) << std::endl; \
  }									\

// fill in the function to evaluate polylogs
#define dilog(x) gsl_sf_dilog(x)


// argument list for 2->2 amplitudes
#define AMP_2_2_ARGS const PS_2_2& ps
// argument list for 2->3 amplitudes
#define AMP_2_3_ARGS const PS_2_3& ps


// this defines references to access Phi_1
#define AMP_H_REFS				\
  double const& fH = FH_eff_1;			\
  double const& fA = FA_eff_1;			\
  double const& At = At_1;			\
  double const& Bt = Bt_1;			\
  double const& mH = M_1;			\
  double const& GammaH = G_1;			\
  double const& mH2 = M2_1;			\
  double const& GammaH2 = G2_1;			\
  
#define AMP_H1_REFS				\
  double const& fH1 = FH_eff_1;			\
  double const& fA1 = FA_eff_1;			\
  double const& At1 = At_1;			\
  double const& Bt1 = Bt_1;			\
  double const& mH1 = M_1;			\
  double const& GammaH1 = G_1;			\

#define AMP_H2_REFS				\
  double const& fH2 = FH_eff_2;			\
  double const& fA2 = FA_eff_2;			\
  double const& At2 = At_2;			\
  double const& Bt2 = Bt_2;			\
  double const& mH2 = M_2;			\
  double const& GammaH2 = G_2;			\
  
#define AMP_H12_REFS				\
  AMP_H1_REFS					\
  AMP_H2_REFS					\

// this defines pointers to the 4 vectors in PS
// here it is assumed that p2 is eliminated from the amplitude
#define AMP_2_2_4VEC_REFS(PS)			\
  using namespace Constants;			\
  using namespace RunParameters;		\
  using namespace Prefactors;			\
  using namespace Bosons;			\
  double AlphaS3 = AlphaS2*AlphaS;		\
  double const& At = At_1;			\
  double const& Bt = Bt_1;			\
  double const& FH0 = FH_1;			\
  double const& FA0 = FA_1;			\
  double const& mH  = M_1;			\
  double const& GammaH = G_1;			\
  double const& LN_MUR2 = LNMU2;		\
  double PREF_V_CA_At = PREF_V_CA*At;		\
  double PREF_V_CA_Bt = PREF_V_CA*Bt;		\
  double PREF_V_CF_At = PREF_V_CF*At;		\
  double PREF_V_CF_Bt = PREF_V_CF*Bt;		\
  double PREF_V_CFCA2_At = PREF_V_CFCA2*At_1;	\
  double PREF_V_CFCA2_Bt = PREF_V_CFCA2*Bt_1;	\
  double PREF_B_At = PREF_B_PHIxQCD*At;		\
  double PREF_B_Bt = PREF_B_PHIxQCD*Bt;		\
  const FV& p1 = PS.p1();			\
  const FV& p2 = PS.p2();			\
  const FV& k1 = PS.k1();			\
  const FV& k2 = PS.k2();			\

#define AMP_2_3_4VEC_REFS(PS)			\
  using namespace Constants;			\
  using namespace RunParameters;		\
  using namespace Prefactors;			\
  using namespace Bosons;			\
  double PREF_R_CA_At = PREF_R_CA*At_1;		\
  double PREF_R_CA_Bt = PREF_R_CA*Bt_1;		\
  double PREF_R_CF_At = PREF_R_CF*At_1;		\
  double PREF_R_CF_Bt = PREF_R_CF*Bt_1;		\
  double PREF_R_CFCA2_At = PREF_R_CFCA2*At_1;	\
  double PREF_R_CFCA2_Bt = PREF_R_CFCA2*Bt_1;	\
  double const& At = At_1;			\
  double const& Bt = Bt_1;			\
  double const& FH0 = FH_1;			\
  double const& FA0 = FA_1;			\
  double const& mH  = M_1;			\
  double const& GammaH = G_1;			\
  double AlphaS3 = AlphaS2*AlphaS;		\
  FV const& p1 = PS.p1();			\
  FV const& p2 = PS.p2();			\
  FV const& k1 = PS.k1();			\
  FV const& k2 = PS.k2();			\
  FV const& p3 = PS.p3();			\

#endif
