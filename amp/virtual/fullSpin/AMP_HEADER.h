

#ifndef BORN_AMP_HEADER
#define BORN_AMP_HEADER

#define AMP_ARGS PS_2_2 const& ps

#define AMP_DEFINITIONS				\
  using namespace Constants;			\
  using namespace Prefactors;			\
  using namespace Bosons;			\
  using namespace RunParameters;		\
  FV const& p1 = ps.p1();			\
  FV const& p2 = ps.p2();			\
  FV const& k1 = ps.k1();			\
  FV const& k2 = ps.k2();			\
  FV const& s1 = ps.s1();			\
  FV const& s2 = ps.s2();			\
  double const& s      = ps.get_s();		\
  double const& y      = ps.get_y();		\
  double const& beta   = ps.get_beta();		\
  double const& beta_y = ps.get_beta_y();	\



#endif
