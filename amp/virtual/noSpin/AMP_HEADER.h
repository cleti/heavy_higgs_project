

#ifndef BORN_AMP_HEADER
#define BORN_AMP_HEADER

#define AMP_ARGS				\
  PS_2_2 const& ps,				\
    AmplitudePrefactors const& ap,		\
    HiggsPrefactors const& hp			\
    

#define AMP_DEFINITIONS				\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  double const& s      = ps.get_s();		\
  double const& y      = ps.get_y();		\
  double const& beta   = ps.get_beta();		\
  double const& beta_y = ps.get_beta_y();	\




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




#define AP_REFS_B(AP)							\
  double const& PREF_B_PHIxPHI      =     AP.PREF_B_PHIxPHI;		\
  double const& PREF_B_PHIxQCD      =     AP.PREF_B_PHIxQCD;		\
  double const& PREF_B_QCDxQCD      =     AP.PREF_B_QCDxQCD;		\
  double const& PREF_B_QCDxQCD_CF   =     AP.PREF_B_QCDxQCD_CF;		\
  double const& PREF_B_QCDxQCD_CA   =     AP.PREF_B_QCDxQCD_CA;		\
  double const& PREF_B_QCDxQCD_CFCA =     AP.PREF_B_QCDxQCD_CFCA;	\



#define AP_REFS_V(AP)					\
  double const& PREF_V_CF    = AP.PREF_V_PHIxQCD_CF;	\
  double const& PREF_V_CA    = AP.PREF_V_PHIxQCD_CA;	\
  double const& PREF_V_CFCA2 = AP.PREF_V_PHIxQCD_CFCA2;	\
  double const& PREF_V_Nf    = AP.PREF_V_PHIxQCD_Nf;	\
  double const& PREF_V_CT    = AP.PREF_V_PHIxQCD_CT;	\
  double const& PREF_V       = AP.PREF_V_PHIxQCD;	\
  double const& PREF_V_PHI    = AP.PREF_V_PHIxPHI;	\
  double const& PREF_V_PHI_CA = AP.PREF_V_PHIxPHI_CA;	\
  double const& PREF_V_PHI_CF = AP.PREF_V_PHIxPHI_CF;	\



#define AP_REFS_R(AP)					\
  double const& PREF_R       = AP.PREF_R_PHIxQCD;	\
  double const& PREF_R_CF    = AP.PREF_R_PHIxQCD_CF;	\
  double const& PREF_R_CA    = AP.PREF_R_PHIxQCD_CA;	\
  double const& PREF_R_CFCA2 = AP.PREF_R_PHIxQCD_CFCA2;	\
  double const& PREF_R_PHI    = AP.PREF_R_PHIxPHI;	\
  double const& PREF_R_PHI_CA = AP.PREF_R_PHIxPHI_CA;	\
  double const& PREF_R_PHI_CF = AP.PREF_R_PHIxPHI_CF;	\


#endif
