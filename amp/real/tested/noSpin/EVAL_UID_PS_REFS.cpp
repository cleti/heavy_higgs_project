
  // normal 2->3 phase space
  FV const& k1 = ps.k1();
  FV const& k2 = ps.k2();
  FV const& p1 = ps.p1();
  FV const& p2 = ps.p2();
  FV const& p3 = ps.p3();

  // reduced dipole phase space
  FV const& P1 = ps_red.p1();
  FV const& P2 = ps_red.p2();
  FV const& K1 = ps_red.k1();
  FV const& K2 = ps_red.k2();

  using namespace Constants;			
  using namespace RunParameters;		
  using namespace Prefactors;			
  using namespace Bosons;			
  double AlphaS3 = AlphaS2*AlphaS;		
  double const& At = At_1;			
  double const& Bt = Bt_1;			
  double const& FH0 = FH_1;			
  double const& FA0 = FA_1;			
  double const& mH  = M_1;			
  double const& GammaH = G_1;					
  double PREF_R_CA_At = PREF_R_CA*At;		
  double PREF_R_CA_Bt = PREF_R_CA*Bt;		
  double PREF_R_CF_At = PREF_R_CF*At;		
  double PREF_R_CF_Bt = PREF_R_CF*Bt;		
  double PREF_R_CFCA2_At = PREF_R_CFCA2*At_1;	
  double PREF_R_CFCA2_Bt = PREF_R_CFCA2*Bt_1;		
