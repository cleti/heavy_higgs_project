

#ifndef REAL_AMP_HEADER
#define REAL_AMP_HEADER

#define AMP_ARGS				\
  const PS_2_3& ps,				\
  const AmplitudePrefactors& ap,		\
  const HiggsPrefactors& hp			\


#define AMP_DEFINITIONS				\
  using namespace Constants;			\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  FV const& p3 = ps.p3();			\




#define HP_REFS_PHIxQCD(HP)			\
  double const& At_fH_re = HP.At_fH_re;		\
  double const& At_fA_re = HP.At_fA_re;		\
  double const& Bt_fH_re = HP.Bt_fH_re;		\
  double const& Bt_fA_re = HP.Bt_fA_re;		\
  double const& At_fH_im = HP.At_fH_im;		\
  double const& At_fA_im = HP.At_fA_im;		\
  double const& Bt_fH_im = HP.Bt_fH_im;		\
  double const& Bt_fA_im = HP.Bt_fA_im;		\



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

#endif
