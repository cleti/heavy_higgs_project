
#include "../inc/PhaseSpace.h"






////// class PS_2   ///////////////////////////////////////////////////////////////////////////////////////
// set incoming particle momenta
void PS_2::set_initial_state(double const& rs)
{
#ifdef DEBUG
  if ( std::isnan(rs) || std::isinf(rs) )
    {
      WARNING("rs is inf/NaN.");
      return;
    }
#endif
      double E  = 0.5*rs;
      p[0][0] =  E;
      p[0][3] =  E;
  
      p[1][0] =  E;
      p[1][3] = -E;
}
void PS_2::swap()
{
  std::swap(p[0],p[1]);
  std::swap(P[0],P[1]);
}
void PS_2::print() const
{
  std::cout << std::endl << "--------------------------------------------------------";
  std::cout << std::endl << "\"" << d_name << " " << this << "\" <--- " << d_parent;
  std::cout << std::endl << "--------------------------------------------------------";
  std::cout << std::endl << " proton momenta [in parton frame]";
  PRINT_4VEC(P[0]);
  PRINT_4VEC(P[1]);
  std::cout << std::endl << "--------------------------------------------------------";
  std::cout << std::endl << "[IS]";
  std::cout << std::endl << "--------------------------------------------------------";
  PRINT_4VEC(p[0]);std::cout << " [ m^2 = " << MSQ(p[0]) << "]" << std::endl;
  PRINT_4VEC(p[1]);std::cout << " [ m^2 = " << MSQ(p[1]) << "]" << std::endl;
}
PS_2::~PS_2()
{
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_1 ///////////////////////////////////////////////////////////////////////////////////////
PS_2_1::PS_2_1(const std::string& nm) : 
  PS_2(nm),
  k(0.0),
  S(0.0),
  d_msq(1.0),
  d_child(nullptr),
  d_x   (0.0)
{}
PS_2_1::PS_2_1(double const& msq, const std::string& nm) : 
  PS_2(nm),
  k(0.0),
  S(0.0),
  d_msq(msq),
  d_child(nullptr),
  d_x   (0.0)
{}
PS_2_1::~PS_2_1()
{
  // remove possible references to this instance
  if (d_child != nullptr) d_child->set_parent(nullptr);
}
int PS_2_1::set()
{
  static double msq_t = -1.0;
  if (d_msq != msq_t)
    {
      set_initial_state(sqrt(d_msq));
      k = p[0] + p[1];
    }
  d_wgt = Constants::TwoPi;
  return 1;
}
void PS_2_1::print() const
{
  PS_2::print();
  std::cout << std::endl << "[FS] -- m1^2 = " << d_msq;
  std::cout << std::endl << "--------------------------------------------------------";
  PRINT_4VEC(k);
  PRINT_4VEC(S);
}
void PS_2_1::FillDistributions(
			       std::vector<HistArray*> & dist,
			       int id,
			       double const & wgt) const
{
  std::cout << "\n FillDistributions: no distributions in case PS 2->1 !!! \n"; exit(1);
}
////// class PS_2_2 ///////////////////////////////////////////////////////////////////////////////////////
PS_2_2::PS_2_2(const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0},
  S{0.0,0.0},
  S_r{0.0,0.0},
  d_msq{1.0,1.0},
  d_beta{0.0,0.0},
  d_child{nullptr,nullptr},
  d_y   (0.0),
  d_phi (0.0),
  d_t11 (0.0),
  d_t12 (0.0),
  d_x   (0.0),
  d_beta_y(0.0)
{}
PS_2_2::PS_2_2(double const& m1sq, double const& m2sq,const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0},
  S{0.0,0.0},
  S_r{0.0,0.0},
  d_msq{m1sq,m2sq},
  d_beta{0.0,0.0},
  d_child{nullptr,nullptr},
  d_y   (0.0),
  d_phi (0.0),
  d_t11 (0.0),
  d_t12 (0.0),
  d_x   (0.0),
  d_beta_y(0.0)
{}
PS_2_2::~PS_2_2()
{
  // remove possible references to this instance
  if (d_child[0] != nullptr) d_child[0]->set_parent(nullptr);
  if (d_child[1] != nullptr) d_child[1]->set_parent(nullptr);
}

int PS_2_2::boost_initial_state()
{
  return boost_initial_state(d_x);
}
int PS_2_2::boost_initial_state(double const& x)
{
  if (x<=0.0) return 0;

  double const& rs   = d_rs;
  double const& rs_x = sqrt(x)*rs;

  return set(rs_x, d_y, d_phi);
}

// boost final state 4-vectors to c.m.f. of the parent phase space
// it is assumed that the 4-vectors are defined in the rest-frame of the parent 4-vector that points to this phase space (the decaying 4-vector)
int PS_2_2::boost_final_state()
{
  if (d_parent != nullptr)
    {
      static LT boost;
      int i = 0;
      while (RANGE(i,2) && d_parent->get_child(i)!=this) ++i;// get index of the 4-vector in parent phase space
      FV const& K = d_parent->get_k(i);

      // boost from z.m.f. into K-restframe
      if (!boost.set_boost(K))
	{
	  WARNING("could not compute boosted 1->2 final state vectors.");
	  SLEEP(3);
	  return 0;
	}

      // boost.apply(S[0]);
      // boost.apply(S[1]);

      // boost from K-restframe into z.m.f.
      boost.invert();
      boost.apply(k[0]);
      boost.apply(k[1]);

      return 1;
    }
  else
    {
      return 0;
    }
}

// set the 4-vector components corresponding to p1,p2 -> k1,k2
// evaluated in the p1-p2/k1-k2 c.m.f., where p1 defines the z-direction !!!
// independent phase space variables are rs=sqrt((p1+p2)^2)
// and the scattering angles y=cos(theta),phi
int PS_2_2::set(double const& rs,
		double const& y,
		double const& phi)
{
  // check if the setting is physical sqrt(s)>m1+m2
  if ( (rs-sqrt(d_msq[0])-sqrt(d_msq[1]))<EPS_BETA2) return 0;
  //////////////////////////////////////////////////////////////
  PS_2::set_initial_state(rs);
  //////////////////////////////////////////////////////////////

  double dm2= d_msq[0]-d_msq[1];
  double P  = 0.5*sqrt(lambda(rs*rs,d_msq[0],d_msq[1]))/rs;
  double E1 = rs - 0.5*(rs-dm2/rs);
  double E2 = 0.5*(rs-dm2/rs);
  double sy = sqrt(1.0-y*y);

#ifdef DEBUG
  if ( std::isnan(P) || std::isinf(P) )
    {
      WARNING("P is inf/NaN, this should better not happen too often!");
      PRINT(lambda(rs*rs,d_msq[0],d_msq[1]));
      SLEEP(3);
      return 0;
    }
  if ( std::isnan(E1) || std::isinf(E1) )
    {
      WARNING("E1 is inf/NaN, too small value for sqrt(s)?");
      PRINT(rs);
      SLEEP(3);
      return 0;
    }
  if ( std::isnan(sy) || std::isinf(sy) )
    {
      WARNING("sy is inf/NaN, incorrect value for y!");
      PRINT(1.0-y*y);
      SLEEP(3);
      return 0;
    }
#endif
  
  // set class members
  double cphi = cos(phi);
  double sphi = sqrt(1.0-cphi*cphi);
  if (phi>Constants::Pi) sphi *= (-1);

  
  k[0][0] = E1;
  k[0][1] = P*sy*cphi;
  k[0][2] = P*sy*sphi;
  k[0][3] = P*y;

  k[1][0] = E2;
  k[1][1] = -P*sy*cphi;
  k[1][2] = -P*sy*sphi;
  k[1][3] = -P*y;


  // copy new values to class member variables
  d_rs   = rs;
  d_beta[0]= P/E1;
  d_beta[1]= P/E2;
  d_y    = y;
  d_phi  = phi;
  d_t11  = 2.0*sp(p[0],k[0]);
  d_t12  = 2.0*sp(p[0],k[1]);
  d_beta_y = d_y*d_beta[0];

  // store store phase space density/weight
  d_wgt = P/(4.0*rs*pow(Constants::TwoPi,2));

  return 1;
}

void PS_2_2::set()
{
  // this function assumes that the vectors p1,p2,k1,k2 have been set up correctly
  // no checks at all !!!
  double s  = (2.0*sp(p[0],p[1]));
  d_rs      = sqrt(s);

  // use invariant description to define beta
  // correpsonds to P/E in the c.m.f. of k1 and k2 !!!
  d_beta[0]= sqrt(1.0-4.0*d_msq[0]/s);
  d_beta[1]= sqrt(1.0-4.0*d_msq[1]/s);
  d_t11  = 2.0*sp(p[0],k[0]);
  d_t12  = 2.0*sp(p[0],k[1]);
  // sp(k1,p2)-sp(k1,p1) = sp(k2,p1) - sp(k1,p1)  = s/2*beta*y
  d_beta_y = (d_t12-d_t11)/s;
  d_y      = d_beta_y/d_beta[0];
  // not needed at the moment ///
  d_phi    = 0.0;
  ///////////////////////////////
}


double PS_2_2::cmp_wgt()
{
  double m1sq = sp(k[0],k[0]);
  double m2sq = sp(k[1],k[1]);
  double s = 2.0*sp(p[0],p[1]);
  return 0.5*sqrt(lambda(s,m1sq,m2sq))/(4.0*s*pow(Constants::TwoPi,2));
}




void PS_2_2::print() const
{
  PS_2::print();
  std::cout << std::endl << "[FS] -- m1^2 = " << d_msq[0] << ", m2^2 = " << d_msq[1];
  std::cout << std::endl << "--------------------------------------------------------";
  PRINT_4VEC(k[0]);
  PRINT_4VEC(k[1]);
  PRINT_4VEC(S[0]);
  PRINT_4VEC(S[1]);
}
void PS_2_2::FillDistributions(
			       std::vector<HistArray*> & dist,
			       int id,
			       double const & wgt) const
{
  // need the lab frame 4-vectors K1,K2 of top/antitop
  static LT boost_to_lab_frame;
  // boost to the rest frame of the protons :
  // P1+P2 = 1/x1*p1 + 1/x2*p2 = lab frame

  boost_to_lab_frame.set_boost(P[0]+P[1],0);
  FV K1 = k[0];
  FV K2 = k[1];
  boost_to_lab_frame.apply(K1);
  boost_to_lab_frame.apply(K2);

  // if (id==H_NLO_R)
  //   {
  //     PRINT_4VEC(P[0]);
  //     PRINT_4VEC(P[1]);
  //     PRINT_4VEC(K1);
  //     PRINT_4VEC(K2);
  //     exit(1);
  //   }
  // tt distributions
  dist[0]->FillOne(id,obs_M12(K1,K2),wgt);
  dist[1]->FillOne(id,obs_PT(K1)    ,wgt);
  // dist[2]->FillOne(id,obs_PT(K2)    ,wgt);
  // need lab frame vectors for these
  dist[3]->FillOne(id,obs_Y(K1)     ,wgt);
  // dist[4]->FillOne(id,obs_Y(K2)     ,wgt);
  dist[5]->FillOne(id,obs_DY(K1,K2) ,wgt);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_3 ///////////////////////////////////////////////////////////////////////////////////////
PS_2_3::PS_2_3(const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0,0.0},
  S{0.0,0.0,0.0},
  S_r{0.0,0.0,0.0},
  d_msq{1.0,1.0,0.0},
  d_beta{0.0,0.0,0.0}, 
  d_child{nullptr,nullptr,nullptr},
  d_y_cm (0.0),
  d_M12  (0.0),
  d_y_12 (0.0),
  d_phi_12(0.0)
{}
PS_2_3::PS_2_3(double const& m1sq, double const& m2sq, double const& m3sq, const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0,0.0},
  S{0.0,0.0,0.0},
  S_r{0.0,0.0,0.0},
  d_msq{m1sq,m2sq,m3sq},
  d_beta{0.0,0.0,0.0}, 
  d_child{nullptr,nullptr,nullptr},
  d_y_cm (0.0),
  d_M12  (0.0),
  d_y_12 (0.0),
  d_phi_12(0.0)
{}
PS_2_3::~PS_2_3()
{
  // remove possible references to this instance
  if (d_child[0] != nullptr) d_child[0]->set_parent(nullptr);
  if (d_child[1] != nullptr) d_child[1]->set_parent(nullptr);
}

// boost final state 4 vectors to c.m.f. of the parent phase space
int PS_2_3::boost_to_parent()
{
  if (d_parent == nullptr) return 0;

  static LT boost;
  int i = 0;
  while (RANGE(i,3) && d_parent->get_child(i)!=this) ++i;
  FV const& K = d_parent->get_k(i);

  if (!boost.set_boost(K))
    {
      WARNING("could not compute boosted 1->2 final state vectors.");
      SLEEP(3);
      return 0;
    }
  boost.invert();
  boost.apply(k[0]);
  boost.apply(k[1]);
  boost.apply(k[2]);
  boost.apply(S[0]);
  boost.apply(S[1]);
  boost.apply(S[2]);
  return 1;
}

int  PS_2_3::set(double const& rs,
		 double const& y_cm,
		 double const& phi_cm,
		 double const& M12,
		 double const& y_12,
		 double const& phi_12)
{
  /////////////////////// split p1 p2 -> k1 k2 k3 into two 2->2 phase spaces
  static PS_2_2 ps_Q3(0.0     ,d_msq[2],"p1 p2 -> Q  k3 [static in PS_2_3::set]"); 
  static PS_2_2 ps_12(d_msq[0],d_msq[1],"Q     -> k1 k2 [static in PS_2_3::set]"); // 
  static LT boost;

  // the first mass serves as integration variable, different in every call
  ps_Q3.set_msq(pow(M12,2),d_msq[2]);
  // set p1,p2 -> Q ,k3 subspace (can be rotated such that phi=0)
  if (!ps_Q3.set(rs,y_cm,phi_cm))
    {
      return 0;
    }

  // set Q -> k1 , k2 subspace
  if (!ps_12.set(M12,y_12,phi_12))
    {
      return 0;
    }

  // boost from Q restframe to current frame
  if (!boost.set_boost(ps_Q3.k1(),1))
    {
      WARNING("could not compute boosted 1->2 final state vectors.");
      SLEEP(3);
      return 0;
    }

  boost.apply(ps_12.k1());
  boost.apply(ps_12.k2());

  // move 4 vectors
  // .. initial state
  p[0] = std::move(ps_Q3.p1());
  p[1] = std::move(ps_Q3.p2());
  // .. final state
  k[0] = std::move(ps_12.k1());
  k[1] = std::move(ps_12.k2());
  k[2] = std::move(ps_Q3.k2());

  // copy new values to class member variables
  d_rs    = rs;
  // d_beta[0] = P1/E1;
  // d_beta[1] = P2/E2;
  // d_beta[2] = P3/E3;
  d_y_cm  = y_cm;
  d_M12   = M12;
  d_y_12  = y_12;
  d_phi_12= phi_12;

  // store phase space density/weight
  d_wgt = 2.0*M12/Constants::TwoPi*ps_12.get_wgt()*ps_Q3.get_wgt();

  return 1;
}



void PS_2_3::print() const
{
  PS_2::print();
  std::cout << std::endl << "[FS]";
  std::cout << std::endl << "--------------------------------------------------------";
  PRINT_4VEC(k[0]);
  PRINT_4VEC(k[1]);
  PRINT_4VEC(k[2]);
  PRINT_4VEC(S[0]);
  PRINT_4VEC(S[1]);
  PRINT_4VEC(S[2]);
}
void PS_2_3::FillDistributions(
			       std::vector<HistArray*> & dist,
			       int id,
			       double const & wgt) const
{
  // need the lab frame 4-vectors K1,K2 of top/antitop
  static LT boost_to_lab_frame;
  // boost to the rest frame of the protons :
  // P1+P2 = 1/x1*p1 + 1/x2*p2 = lab frame
  boost_to_lab_frame.set_boost(P[0]+P[1],0);
  FV K1 = k[0];
  FV K2 = k[1];
  boost_to_lab_frame.apply(K1);
  boost_to_lab_frame.apply(K2);
  
  // tt distributions
  dist[0]->FillOne(id,obs_M12(K1,K2),wgt);
  dist[1]->FillOne(id,obs_PT(K1)    ,wgt);
  // dist[2]->FillOne(id,obs_PT(K2)    ,wgt);
  // need lab frame vectors for these
  dist[3]->FillOne(id,obs_Y(K1)     ,wgt);
  // dist[4]->FillOne(id,obs_Y(K2)     ,wgt);
  dist[5]->FillOne(id,obs_DY(K1,K2) ,wgt);
}
