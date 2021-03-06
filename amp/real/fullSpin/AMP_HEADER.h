

#ifndef REAL_AMP_HEADER
#define REAL_AMP_HEADER

#define AMP_ARGS				\
  PS_2_3 const& ps,				\
    AmplitudePrefactors const& ap,		\
    HiggsPrefactors const& hp			\

#define AMP_DEFINITIONS				\
  using namespace Constants;			\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  FV const& p3 = ps.p3();			\
  FV const& s1 = ps.s1();			\
  FV const& s2 = ps.s2();			\



#define HP_REFS_PHIxPHI(HP)				\
  double const& At2_fH2_De = HP.At2_fH2_De;		\
  double const& At2_fA2_De = HP.At2_fA2_De;		\
  double const& Bt2_fH2_De = HP.Bt2_fH2_De;		\
  double const& Bt2_fA2_De = HP.Bt2_fA2_De;		\
  double const& At_Bt_fH2_De = HP.At_Bt_fH2_De;		\
  double const& At_Bt_fA2_De = HP.At_Bt_fA2_De;		\
  double const& At_Bt_fH2_DeIM = HP.At_Bt_fH2_DeIM;	\
  double const& At_Bt_fA2_DeIM = hp.At_Bt_fA2_DeIM;	\


#define AP_REFS_R(AP)					\
  double const& PREF_R       = AP.PREF_R_PHIxQCD;	\
  double const& PREF_R_CF    = AP.PREF_R_PHIxQCD_CF;	\
  double const& PREF_R_CA    = AP.PREF_R_PHIxQCD_CA;	\
  double const& PREF_R_CFCA2 = AP.PREF_R_PHIxQCD_CFCA2;	\
  double const& PREF_R_PHI    = AP.PREF_R_PHIxPHI;	\
  double const& PREF_R_PHI_CA = AP.PREF_R_PHIxPHI_CA;	\
  double const& PREF_R_PHI_CF = AP.PREF_R_PHIxPHI_CF;	\
  double const& PREF_V_CA    = AP.PREF_V_PHIxQCD_CA;	\
  
#endif
