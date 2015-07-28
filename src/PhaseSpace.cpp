
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
  P[0].swap(P[1]);
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
  
#ifdef WITH_T_SPIN
  ///////////////////////////////////////////////////////
  // spin-dependent observables /////////////////////////
  ///////////////////////////////////////////////////////
  for (auto dist: dist_vec)
    {
      HistArray* hist = dist->GetHistograms();
      double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
      double obs_val  = (*dist)(this)*mass_dim;
      double avg_val  = dist->Avg(this);
      dist->GetHistograms()->FillOne(id,obs_val,avg_val*wgt);
    }
  ///////////////////////////////////////////////////////
#else
  ///////////////////////////////////////////////////////
  // spin-independent observables ///////////////////////
  ///////////////////////////////////////////////////////  
  // for observables like rapidity we need the lab frame 4-vectors K1,K2 of top/antitop
  static LT boost_to_lab_frame;
  boost_to_lab_frame.set_boost(P[0]+P[1],0);
  // make copies of the top/antitop 4-vectors before the boost,
  // the vectors of this ps might be used somewhere else!!!
  static PS_2_2 ps_lab("ps 2->2, lab frame, static in PS_2_2::FillDistributions");
  ps_lab = *this;
  boost_to_lab_frame.apply(ps_lab.k1());
  boost_to_lab_frame.apply(ps_lab.k2());

  
  for (auto dist: dist_vec)
    {
      HistArray* hist = dist->GetHistograms();
      double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
      double obs_val  = (*dist)(&ps_lab)*mass_dim;
      double avg_val  = dist->Avg(&ps_lab);
      dist->GetHistograms()->FillOne(id,obs_val,avg_val*wgt);
      //e.first->FillOne(id,e.second(&ps_lab)*m,wgt);
    }
  ///////////////////////////////////////////////////////
#endif
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
#ifdef WITH_T_SPIN
  ///////////////////////////////////////////////////////
  // spin-dependent observables /////////////////////////
  ///////////////////////////////////////////////////////  
  // to be implemented ...
  ERROR("the code for distributions on PS_2_3 for polarized amplitudes is still missing!");
  ///////////////////////////////////////////////////////
#else
  ///////////////////////////////////////////////////////
  // spin-independent observables ///////////////////////
  ///////////////////////////////////////////////////////  
  // for observables like rapidity we need the lab frame 4-vectors K1,K2 of top/antitop
  static LT boost_to_lab_frame;
  boost_to_lab_frame.set_boost(P[0]+P[1],0);
  // // make copies of the top/antitop 4-vectors before the boost,
  // // the vectors of this ps might be used somewhere else!!!
  // FV K1 = k[0];
  // FV K2 = k[1];
  // boost_to_lab_frame.apply(K1);
  // boost_to_lab_frame.apply(K2);
  // make copies of the top/antitop 4-vectors before the boost,
  // the vectors of this ps might be used somewhere else!!!
  static PS_2_3 ps_lab("ps 2->3, lab frame, static in PS_2_3::FillDistributions");
  ps_lab = *this;
  boost_to_lab_frame.apply(ps_lab.k1());
  boost_to_lab_frame.apply(ps_lab.k2());
  

	  
  
  // // if (obs_PT12(K1,K2)*mScale < 1.0)
  // //   {
  // //     PRINT(obs_PT12(K1,K2)*mScale);
  // //     PRINT(k[2][0]*mScale); // gluon energy
  // //     PRINT((p[0][0]*k[2][0]-sp(p[0],k[2]))/(p[0][0]*k[2][0]));  // gluon scattering angle
  // //   }
  
  // dist[0].first->FillOne(id,obs_M12(K1,K2)*mScale,wgt);
  // dist[1].first->FillOne(id,obs_PT(K1)*mScale    ,wgt);
  // // dist[2]->FillOne(id,obs_PT(K2)    ,wgt);
  // dist[2].first->FillOne(id,obs_PT12(K1,K2)*mScale,wgt);
  // // need lab frame vectors for these
  // dist[3].first->FillOne(id,obs_Y(K1),wgt);
  // dist[4].first->FillOne(id,obs_Y(K2),wgt);
  // // dist[6].first->FillOne(id,obs_DY(K1,K2) ,wgt);



  for (auto dist: dist_vec)
    {
      HistArray* hist = dist->GetHistograms();
      double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
      double obs_val  = (*dist)(&ps_lab)*mass_dim;
      double avg_val  = dist->Avg(&ps_lab);
      dist->GetHistograms()->FillOne(id,obs_val,avg_val*wgt);
    }
  
  // for (auto dist: dist_vec)
  //   {
  //     HistArray* hist = dist->GetHistograms();
  //     double mass_dim = (hist->GetMassDim())?std::pow(mScale,hist->GetMassDim()):1.0;
  //     double obs_val  = (*dist)(&ps_lab)*mass_dim;
  //     double avg_val  = dist->Avg(&ps_lab);
  //     dist->GetHistograms()->FillOne(id,obs_val,avg_val);
  //     //e.first->FillOne(id,e.second(&ps_lab)*m,wgt);
  //     // PRINT(mass_dim);
  //     // PRINT(obs_val);
  //     // PRINT(avg_val);
  //   }
  //EXIT(1);
  ///////////////////////////////////////////////////////
#endif
}





/////////////////////////////////////////////////////////////////////////////////////////////
// definition of observables used for this project, feel free to add more ///////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// mass of 2-particle system [GeV]
double obs_M12(FV const& k1, FV const& k2)
{
  return sqrt(sp(k1,k1)+sp(k2,k2)+2.0*sp(k1,k2));
}
double OBS_M12(const PS_2* ps)
{
  return obs_M12(ps->k1(),ps->k2());
}


// transverse momentum (transverse plane is the x-y-plane) 
double obs_PT(FV const& k)
{
  return sqrt(k[1]*k[1]+k[2]*k[2]);
}

double OBS_PT1(const PS_2* ps)
{
  return obs_PT(ps->k1());
}
double OBS_PT2(const PS_2* ps)
{
  return obs_PT(ps->k2());
}



// transverse momentum of two particle system |k_{T,1}+k_{T,2}| 
double obs_PT12(FV const& k1, FV const& k2)
{
  return sqrt(pow(k1[1]+k2[1],2)+pow(k1[2]+k2[2],2)); 
}
double OBS_PT12(const PS_2* ps)
{
  //PRINT(obs_PT12(ps->k1(),ps->k2()));
  return obs_PT12(ps->k1(),ps->k2());
}



// pseudo-rapidity
double obs_Y(FV const& k)
{
  return 0.5*log(((k[0]+k[3])/(k[0]-k[3])));
}
double OBS_Y1(const PS_2* ps)
{
  //PRINT(obs_Y(ps->k1()));
  return obs_Y(ps->k1());
}
double OBS_Y2(const PS_2* ps)
{
  //PRINT(obs_Y(ps->k2()));
  return obs_Y(ps->k2());
}



// difference of abs. values of transverse rapidities
double obs_DY(FV const& k1, FV const& k2)
{
  return fabs(obs_Y(k1))-fabs(obs_Y(k2));
}
double OBS_DY12(const PS_2* ps)
{
  // PRINT(obs_DY(ps->k1(),ps->k2()));
  // PRINT_4VEC(ps->k1());
  // PRINT_4VEC(ps->k2());
  return obs_DY(ps->k1(),ps->k2());
}



// 3-dim opening angle of the two vectors
double obs_PHI(FV const& k1, FV const& k2)
{
  // k1 dot k2 / |k1| / |k2|
  double K1 = LEN(k1);
  double K2 = LEN(k2);
  return (k1[1]*k2[1]+k1[2]*k2[2]+k1[3]*k2[3])/(K1*K2);
}
// 2-dim transversal opening angle of the two vectors
double obs_PHIT(FV const& k1, FV const& k2)
{
  // k1_T dot k2_T / |k1_T| / |k2_T|
  double K1 = sqrt(k1[1]*k1[1]+k1[2]*k1[2]);
  double K2 = sqrt(k2[1]*k2[1]+k2[2]*k2[2]);
  return (k1[1]*k2[1]+k1[2]*k2[2])/(K1*K2);
}
// spatial triple product of three 4-vectors k1,k2,k3
double obs_TriProd(FV const& k1, FV const& k2, FV const& k3)
{
  // = (k1 x k2) dot k3 
  double t3 = k1[1] * k2[2] - k1[2] * k2[1];
  double t7 = -k1[1] * k2[3] + k1[3] * k2[1];
  double t11 = k1[2] * k2[3] - k1[3] * k2[2];
  return (t11 * k3[1] + t3 * k3[3] + t7 * k3[2]);
}
double obs_TriProdN(FV const& k1, FV const& k2, FV const& k3)
{
  // = (k1 x k2) dot k3 / |k1 x k2| / |k3|
  double t3 = k1[1] * k2[2] - k1[2] * k2[1];
  double t7 = -k1[1] * k2[3] + k1[3] * k2[1];
  double t11 = k1[2] * k2[3] - k1[3] * k2[2];
  double t14 = t3 * t3;
  double t15 = t7 * t7;
  double t16 = t11 * t11;
  double t18 = std::pow(k3[1], 2);
  double t19 = std::pow(k3[2], 2);
  double t20 = std::pow(k3[3], 2);
  double t23 = std::sqrt((t14 + t15 + t16) * (t18 + t19 + t20));
  return (t11 * k3[1] + t3 * k3[3] + t7 * k3[2]) / t23;
}





#ifdef WITH_T_SPIN  
double OBS_D12(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  return (-3.0)*obs_PHI(ps->s1_r(),ps->s2_r())/Constants::kappa_p/Constants::kappa_m;
}

double OBS_CP1(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1 is the top momentum in the ttbar zero-momentum frame
  return obs_TriProdN(ps->s1_r(),ps->s2_r(),ps->k1())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_CP2(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1 is the top momentum in the ttbar zero-momentum frame
  return obs_TriProdN(ps->s1_r(),ps->s2_r(),ps->k2())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_HEL12(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (-9.0)*obs_PHI(ps->k1(),ps->s1_r())*obs_PHI(ps->k2(),ps->s2_r())/Constants::kappa_p/Constants::kappa_m;
}
double OBS_B1(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (3.0)*obs_PHI(ps->k1(),ps->s1_r())/Constants::kappa_p;
}
double OBS_B2(const PS_2* ps)
{
  // s1_r, s2_r are the lepton momenta in t/tbar rest frames
  // k1, k2 are the top/antitop momenta in the ttbar zero-momentum frame
  return (3.0)*obs_PHI(ps->k2(),ps->s2_r())/Constants::kappa_m;
}
#endif
