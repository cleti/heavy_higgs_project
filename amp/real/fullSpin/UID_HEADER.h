

#ifndef UID_HEADER
#define UID_HEADER

#define UID_ARGS PS_2_3 const& ps,PS_2_2 const& ps_red


#define UID_DEFINITIONS				\
  using namespace Constants;			\
  using namespace Prefactors;			\
  using namespace Bosons;			\
  using namespace RunParameters;		\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  FV const& p3 = ps.p3();			\
  FV const& P1 = ps_red.p1();			\
  FV const& P2 = ps_red.p2();			\
  FV const& K1 = ps_red.k1();			\
  FV const& K2 = ps_red.k2();			\
  FV const& S1 = ps_red.s1();			\
  FV const& S2 = ps_red.s2();			\

#endif
