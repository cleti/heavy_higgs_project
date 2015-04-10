

#ifndef REAL_AMP_HEADER
#define REAL_AMP_HEADER

#define AMP_ARGS PS_2_3 const& ps

#define AMP_DEFINITIONS				\
  using namespace Constants;			\
  using namespace AmpPrefactors;			\
  using namespace HiggsBosons;			\
  using namespace RunParameters;		\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  FV const& p3 = ps.p3();			\
  FV const& s1 = ps.s1();			\
  FV const& s2 = ps.s2();			\

#endif
