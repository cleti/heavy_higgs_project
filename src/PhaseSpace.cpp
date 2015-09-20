
#include "../inc/PhaseSpace.h"






///////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2   ///////////////////////////////////////////////////////////////////////////
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


void PS_2::swap_initial_state()
{
  p[0].swap(p[1]);
  //P[0].swap(P[1]);
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


const FV PS_2::nullvec(0);

///////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_1 ///////////////////////////////////////////////////////////////////////////
PS_2_1::PS_2_1(const std::string& nm) : 
  PS_2(nm),
  k(0.0),
#ifdef WITH_T_SPIN  
  S(0.0),
#endif
  d_msq(1.0),
  d_child(nullptr),
  d_x   (0.0)
{}


PS_2_1::PS_2_1(double const& msq, const std::string& nm) : 
  PS_2(nm),
  k(0.0),
#ifdef WITH_T_SPIN   
  S(0.0),
#endif
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
#ifdef WITH_T_SPIN   
  PRINT_4VEC(S);
#endif
  FV Q = p[0]+p[1]-k[0];
  std::cout << std::endl << "Energy-momentum conservation:";
  PRINT_4VEC(Q);
}


void PS_2_1::FillDistributions(
			       DistVec& dist,
			       H_Index id,
			       double const & wgt,
			       double const& mScale) const
{
  std::cout << "\n FillDistributions: no distributions in case PS 2->1 !!! \n"; exit(1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_2 ///////////////////////////////////////////////////////////////////////////
PS_2_2::PS_2_2(const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0},
#ifdef WITH_T_SPIN  	       
  S{0.0,0.0},
  S_r{0.0,0.0},
#endif	       
  d_msq{1.0,1.0},
  d_beta{0.0,0.0},
  d_child{nullptr,nullptr},
  d_y   (0.0),
  d_phi (0.0),
  d_t11 (0.0),
  d_t12 (0.0),
  d_x   (0.0),
  d_beta_y(0.0)
{
  //  std::cout << std::endl << " PS_2_2 constructor A " << std::endl; 
}
PS_2_2::PS_2_2(double const& m1sq, double const& m2sq,const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0},
#ifdef WITH_T_SPIN  							
  S{0.0,0.0},
  S_r{0.0,0.0},
#endif							
  d_msq{m1sq,m2sq},
  d_beta{0.0,0.0},
  d_child{nullptr,nullptr},
  d_y   (0.0),
  d_phi (0.0),
  d_t11 (0.0),
  d_t12 (0.0),
  d_x   (0.0),
  d_beta_y(0.0)
{
  //  std::cout << std::endl << " PS_2_2 constructor B " << std::endl; 
}
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

// set the 4-vector components for p1,p2 -> k1,k2 process
// as seen in the p1-p2 / k1-k2 z.m.f., where p1 defines the z-direction !!!
// independent phase space variables are the c.m.e. rs=sqrt((p1+p2)^2)
// and the scattering angles y=cos(theta),phi
int PS_2_2::set(double const& rs,
		double const& y,
		double const& phi)
{
  // check if the setting is physical, i.e. sqrt(s)>m1+m2
  if ( (rs-sqrt(d_msq[0])-sqrt(d_msq[1]))<EPS_BETA2) return 0;
  //////////////////////////////////////////////////////////////
  // set up incoming momentum 4-vectors
  PS_2::set_initial_state(rs);
  //////////////////////////////////////////////////////////////

  // compute outgoing particle energies and momentum modulus
  double dm2= d_msq[0]-d_msq[1];
  double P  = 0.5*sqrt(lambda(rs*rs,d_msq[0],d_msq[1]))/rs;
  double E1 = rs - 0.5*(rs-dm2/rs);
  double E2 = 0.5*(rs-dm2/rs);
  double sy = sqrt(1.0-y*y);

#ifdef DEBUG
  if ( std::isnan(P) || std::isinf(P) )
    {
      std::cout << std::endl << get_name() << std::endl;
      WARNING("P is inf/NaN, this should better not happen too often!");
      PRINT(lambda(rs*rs,d_msq[0],d_msq[1]));
      SLEEP(3);
      return 0;
    }
  if ( std::isnan(E1) || std::isinf(E1) )
    {
      std::cout << std::endl << get_name() << std::endl;
      WARNING("E1 is inf/NaN, too small value for sqrt(s)?");
      PRINT(rs);
      SLEEP(3);
      return 0;
    }
  if ( std::isnan(sy) || std::isinf(sy) )
    {
      std::cout << std::endl << get_name() << std::endl;
      WARNING("sy is inf/NaN, incorrect value for y!");
      PRINT(1.0-y*y);
      SLEEP(3);
      return 0;
    }
#endif

  if (phi==0.0)
    {
      // set up the two outgoing momentum 4-vectors
      k[0][0] = E1;
      k[0][1] = P*sy;
      k[0][2] = 0.0;
      k[0][3] = P*y;

      k[1][0] = E2;
      k[1][1] = -P*sy;
      k[1][2] = 0.0;
      k[1][3] = -P*y;
    }
  else
    {
      // polar angle
      double cphi = cos(phi);
      // need to switch sign for phi>pi if we compute sin(phi) in this way
      double sphi = sqrt(1.0-cphi*cphi);
      if (phi>Constants::Pi) sphi *= (-1);

      // set up the two outgoing momentum 4-vectors
      k[0][0] = E1;
      k[0][1] = P*sy*cphi;
      k[0][2] = P*sy*sphi;
      k[0][3] = P*y;

      k[1][0] = E2;
      k[1][1] = -P*sy*cphi;
      k[1][2] = -P*sy*sphi;
      k[1][3] = -P*y;
    }
  
  // store phase space variables and compute invariants
  d_rs   = rs;
  d_beta[0]= P/E1;
  d_beta[1]= P/E2;
  d_y    = y;
  d_phi  = phi;
  d_t11  = 2.0*sp(p[0],k[0]);
  d_t12  = 2.0*sp(p[0],k[1]);
  d_beta_y = d_y*d_beta[0];

  // store store phase space density/weight for later use
  d_wgt = P/(4.0*rs*pow(Constants::TwoPi,2));

  return 1;
}

/*
  this function assumes that the 4-vector components of p1,p2,k1,k2
  have been set up correctly prior to this call and that they are given in
  the p1-p2 / k1-k2 z.m.f. , to get this ps object in a consistent state
  we then have to derive masses, invariants and other phase space variables 
*/
int PS_2_2::set(double const& x)
{
  // no checks at all !!!
  double s  = (2.0*sp(p[0],p[1]));
  
#ifdef DEBUG
  if (s<=0.0)
    {
      WARNING("s is <=0!");
      SLEEP(3);
      return 0;
    }
#endif
  
  d_rs      = sqrt(s);
  d_x       = x;
  double m1sq = sp(k[0],k[0]);
  double m2sq = sp(k[1],k[1]);
  
#ifdef DEBUG
  if (s<=4.0*std::max(m1sq,m2sq))
    {
      std::cout << std::endl << get_name() << std::endl;
      WARNING("s is below threshold given by the final state masses!");
      SLEEP(3);
      return 0;
    }
#endif

  // outgoing particle masses
  d_msq[0] = m1sq;
  d_msq[1] = m2sq;
  // use invariant description to define beta
  // correpsonds to P/E in the c.m.f. of k1 and k2 !!!
  d_beta[0]= sqrt(1.0-4.0*d_msq[0]/s);
  d_beta[1]= sqrt(1.0-4.0*d_msq[1]/s);
  d_t11  = 2.0*sp(p[0],k[0]);
  d_t12  = 2.0*sp(p[0],k[1]);
  // sp(k1,p2)-sp(k1,p1) = sp(k2,p1) - sp(k1,p1)  = s/2*beta*y
  d_beta_y = (d_t12-d_t11)/s;
  d_y      = d_beta_y/d_beta[0];
  // not needed at the moment //
  d_phi    = 0.0;
  return 1;
}

// recompute the phase space weight
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
#ifdef WITH_T_SPIN
  PRINT_4VEC(S[0]);
  PRINT_4VEC(S[1]);
#endif
  FV Q = p[0]+p[1]-k[0]-k[1];
  std::cout << std::endl << "Energy-momentum conservation:";
  PRINT_4VEC(Q);
}


void PS_2_2::FillDistributions(
			       DistVec& dist_vec,
			       H_Index id,
			       double const& wgt,
			       double const& mScale) const
{
  // at the moment we are in the parton z.m.f.
  // for observables like rapidity we need the lab frame 4-vectors
  // in case 2->2: parton z.m.f. = tt z.m.f.
  static LT boost_to_lab_frame;
  
  // lab frame = z.m.f. of P[0] + P[1]
  if (!boost_to_lab_frame.set_boost(P[0]+P[1],0))
    {
      WARNING("could not set boost to lab frame, skipping distributions");
      return;
    }

  // store the lab frame 4-vectors in a new PS_2_3 instance
  static PS_2_2 ps_lab("ps 2->2, lab frame, static in PS_2_2::FillDistributions");
  boost_to_lab_frame.apply_cpy(P[0],ps_lab.P1());
  boost_to_lab_frame.apply_cpy(P[1],ps_lab.P2());
  boost_to_lab_frame.apply_cpy(k[0],ps_lab.k1());
  boost_to_lab_frame.apply_cpy(k[1],ps_lab.k2());
#ifdef WITH_T_SPIN
  boost_to_lab_frame.apply_cpy(S[0],ps_lab.s1());
  boost_to_lab_frame.apply_cpy(S[1],ps_lab.s2());
  // do not transform s1_r / s2_r, they contain the restframe spin 4-vectors.
#endif
  // NOTE: if there are observables depending on others than the above transformed 4-vectors,
  // they must be added here!

  for (auto dist: dist_vec)
    {
      HistArray* hist = dist->GetHistograms();
      double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
      double obs_val  = dist->Obs(&ps_lab,this)*mass_dim;
      double avg_val  = dist->Avg(&ps_lab,this);
      dist->GetHistograms()->FillOne(id,obs_val,avg_val*wgt);
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_3 ///////////////////////////////////////////////////////////////////////////
PS_2_3::PS_2_3(const std::string& nm) : 
  PS_2(nm),
  k{0.0,0.0,0.0},
#ifdef WITH_T_SPIN  									    
  S{0.0,0.0,0.0},
  S_r{0.0,0.0,0.0},
#endif									    
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
#ifdef WITH_T_SPIN  									    
  S{0.0,0.0,0.0},
  S_r{0.0,0.0,0.0},
#endif									    
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
      std::cout << std::endl << get_name() << std::endl;
      WARNING("could not compute boosted 1->2 final state vectors.");
      SLEEP(3);
      return 0;
    }
  boost.invert();
  boost.apply(k[0]);
  boost.apply(k[1]);
  boost.apply(k[2]);
#ifdef WITH_T_SPIN  
  boost.apply(S[0]);
  boost.apply(S[1]);
  boost.apply(S[2]);
#endif 
  return 1;
}

int  PS_2_3::set(
		 double const& rs,
		 double const& y_cm,
		 double const& phi_cm,
		 double const& M12,
		 // //////////////////////////////////////////////
		 // // runs from a12 to 1, used instead of M12
		 // double const& r12,
		 // // lower bound of integrationd variable r12
		 // double const& a12,
		 // //////////////////////////////////////////////
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
      std::cout << std::endl << get_name() << std::endl;
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
#ifdef WITH_T_SPIN    
  PRINT_4VEC(S[0]);
  PRINT_4VEC(S[1]);
  PRINT_4VEC(S[2]);
#endif
  FV Q = p[0]+p[1]-k[0]-k[1]-k[2];
  std::cout << std::endl << "Energy-momentum conservation:";
  PRINT_4VEC(Q);
}

void PS_2_3::FillDistributions(
			       DistVec& dist_vec,
			       H_Index id,
			       double const & wgt,
			       double const& mScale) const
{
  // at the moment we are in the parton z.m.f.
  // for observables like rapidity we need the lab frame 4-vectors
  // as well as the 4-vectors in the tt-z.m.f.
  static LT boost_to_lab_frame;
  static LT boost_to_tt_zmf;
  
  // lab frame = z.m.f. of P[0] + P[1]
  if (!boost_to_lab_frame.set_boost(P[0]+P[1],0))
    {
      WARNING("could not set boost to lab frame, skipping distributions");
      return;
    }

  // tt z.m.f. = z.m.f. of k[0]+k[1]
  if (!boost_to_tt_zmf.set_boost(k[0]+k[1],0))
    {
      WARNING("could not set boost to tt z.m.f., skipping distributions");
      return;
    }

  // FV K12 = k1()+k2();
  // PRINT_4VEC(K12);
  // PRINT_4VEC(k1());
  // PRINT_4VEC(k2());
  // PRINT_4VEC(P1());
  // PRINT_4VEC(P2());
  
  // store the lab frame 4-vectors in a new PS_2_3 instance
  static PS_2_3 ps_lab("ps 2->3, lab frame, static in PS_2_3::FillDistributions");
  boost_to_lab_frame.apply_cpy(P[0],ps_lab.P1());
  boost_to_lab_frame.apply_cpy(P[1],ps_lab.P2());
  boost_to_lab_frame.apply_cpy(k[0],ps_lab.k1());
  boost_to_lab_frame.apply_cpy(k[1],ps_lab.k2());
#ifdef WITH_T_SPIN
  boost_to_lab_frame.apply_cpy(S[0],ps_lab.s1());
  boost_to_lab_frame.apply_cpy(S[1],ps_lab.s2());
#endif

  // K12 = ps_lab.k1()+ps_lab.k2();
  // PRINT_4VEC(K12);
  // PRINT_4VEC(ps_lab.k1());
  // PRINT_4VEC(ps_lab.k2());
  // PRINT_4VEC(ps_lab.P1());
  // PRINT_4VEC(ps_lab.P2());
  
  // store the tt z.m.f. 4-vectors in a new PS_2_3 instance
  static PS_2_3 ps_tt("ps 2->3, tt z.m.f., static in PS_2_3::FillDistributions");  
  // now convert the 4-vectors of this instance to the tt z.m.f.
  boost_to_tt_zmf.apply_cpy(P[0],ps_tt.P1());
  boost_to_tt_zmf.apply_cpy(P[1],ps_tt.P2());
  boost_to_tt_zmf.apply_cpy(k[0],ps_tt.k1());
  boost_to_tt_zmf.apply_cpy(k[1],ps_tt.k2());
#ifdef WITH_T_SPIN
  boost_to_tt_zmf.apply_cpy(S[0],ps_tt.s1());
  boost_to_tt_zmf.apply_cpy(S[1],ps_tt.s2());
#endif  

  // K12 = ps_tt.k1()+ps_tt.k2();
  // PRINT_4VEC(K12);
  // PRINT_4VEC(ps_tt.k1());
  // PRINT_4VEC(ps_tt.k2());
  // PRINT_4VEC(ps_tt.P1());
  // PRINT_4VEC(ps_tt.P2());
  
  // NOTE: only the top and antitop 4-vectors are transformed
  // since we do not consider observables depending on the other 4-vectors

  for (auto dist: dist_vec)
    {
      HistArray* hist = dist->GetHistograms();
      double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
      double obs_val  = dist->Obs(&ps_lab,&ps_tt)*mass_dim;
      double avg_val  = dist->Avg(&ps_lab,&ps_tt);
      dist->GetHistograms()->FillOne(id,obs_val,avg_val*wgt);
    }
}


