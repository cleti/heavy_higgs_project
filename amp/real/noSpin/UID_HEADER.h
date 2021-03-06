

#ifndef UID_HEADER
#define UID_HEADER



#define UID_ARGS				\
  PS_2_3 const& ps,				\
    PS_2_2 const& ps_red,			\
    const HiggsModel& hm



#define UID_Q_QG_ARGS				\
  PS_2_3 const& ps,				\
    PS_2_2 const& ps_red,			\
    const HiggsModel& hm,			\
    FV const& Q,				\
    double const& Vdiag,			\
    double const& Vnd



#define UID_DEFINITIONS						\
  using namespace Constants;					\
  FV const& p1 = ps.p1();					\
  FV const& p2 = ps.p2();					\
  FV const& k1 = ps.k1();					\
  FV const& k2 = ps.k2();					\
  FV const& p3 = ps.p3();					\
  FV const& P1 = ps_red.p1();					\
  FV const& P2 = ps_red.p2();					\
  FV const& K1 = ps_red.k1();					\
  FV const& K2 = ps_red.k2();					\
  const AmplitudePrefactors& ap = hm.GetAmpPrefactors();	\
  const HiggsPrefactors&     hp = hm.GetHiggsPrefactors();	\
  const double& AlphaS2 = hm.AlphaS2();


#ifndef HP_REFS_PHIxQCD
#define HP_REFS_PHIxQCD(HP)			\
  double const& At_fH_re = HP.At_fH_re;		\
  double const& At_fA_re = HP.At_fA_re;		\
  double const& Bt_fH_re = HP.Bt_fH_re;		\
  double const& Bt_fA_re = HP.Bt_fA_re;		\
  double const& At_fH_im = HP.At_fH_im;		\
  double const& At_fA_im = HP.At_fA_im;		\
  double const& Bt_fH_im = HP.Bt_fH_im;		\
  double const& Bt_fA_im = HP.Bt_fA_im;		
#endif


#ifndef HP_REFS_PHIxPHI
#define HP_REFS_PHIxPHI(HP)				\
  double const& At2_fH2_De = HP.At2_fH2_De;		\
  double const& At2_fA2_De = HP.At2_fA2_De;		\
  double const& Bt2_fH2_De = HP.Bt2_fH2_De;		\
  double const& Bt2_fA2_De = HP.Bt2_fA2_De;		\
  double const& At_Bt_fH2_De = HP.At_Bt_fH2_De;		\
  double const& At_Bt_fA2_De = HP.At_Bt_fA2_De;		\
  double const& At_Bt_fH2_DeIM = HP.At_Bt_fH2_DeIM;	\
  double const& At_Bt_fA2_DeIM = HP.At_Bt_fA2_DeIM;
#endif


#ifndef AP_REFS_R
#define AP_REFS_R(AP)					\
  double const& PREF_R       = AP.PREF_R_PHIxQCD;	\
  double const& PREF_R_CF    = AP.PREF_R_PHIxQCD_CF;	\
  double const& PREF_R_CA    = AP.PREF_R_PHIxQCD_CA;	\
  double const& PREF_R_CFCA2 = AP.PREF_R_PHIxQCD_CFCA2;	\
  double const& PREF_R_PHI    = AP.PREF_R_PHIxPHI;	\
  double const& PREF_R_PHI_CA = AP.PREF_R_PHIxPHI_CA;	\
  double const& PREF_R_PHI_CF = AP.PREF_R_PHIxPHI_CF;
#endif



#endif
