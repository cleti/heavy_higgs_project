
/*! \file
  \brief Phase space classes for 2->1, 2->2 and 2->3 scattering processes.
*/ 


#ifndef PHASEPSACE_H
#define PHASEPSACE_H


#include "Global.h"
#include "Makros.h"
#include "Lorentz.h"
#include "HistArray.h"




#define EPS_BETA2 1e-10








inline double lambda (double const& x, double const& y, double const& z)
{
  return x*x+y*y+z*z-2.0*(x*y+x*z+y*z);
}


class PS_Named {
 protected:
  //! name of the PS instance, useful to locate errors
  std::string d_name;
 public:
  PS_Named(std::string const& nm) : d_name(nm) {}
  PS_Named const& operator=(PS_Named const& rhs) { return *this; }
  void set_name(std::string const& name) { d_name = name; }
};

////// class PS_2 [abstract base class] ///////////////////////////////////////////////////////////
/*!
  Abstract base class for 2->X phase space classes, holds among others the 4-vectors of the incoming partons.
  \sa Lorentz.h, HistArray.h
*/
class PS_2: PS_Named {
 protected:
  //! square root of the current partonic c.m.e.
  double d_rs;
  //! incoming parton 4-momenta
  FV p[2];
  //! incoming proton 4-momenta, used to define the boost that connects parton z.m.f. and lab frame
  FV P[2];
  //! current phase space density
  double d_wgt;
  //! not used
  bool   d_decay;
  //! pointer to parent PS instance, used to combine phase spaces 2->i + 2->j = 2->i+j-1
  PS_2*  d_parent;

  void set_initial_state(double const& s);

 public:
  explicit PS_2(const std::string& nm = "") : PS_Named(nm), d_rs(0.0), p{0.0,0.0}, P{0.0,0.0}, d_wgt(0.0), d_decay(false), d_parent(nullptr) {}
  virtual ~PS_2();
  
  //! read/write access to first incoming parton 4-momentum
  FV& p1()  { return p[0]; }
  //! read/write access to second incoming parton 4-momentum
  FV& p2()  { return p[1]; }
  //! read only access to first incoming parton 4-momentum
  FV const& p1() const { return p[0]; }
  //! read only access to second incoming parton 4-momentum
  FV const& p2() const { return p[1]; }
  //! read/write access to first proton 4-momentum
  FV& P1()  { return P[0]; }
  //! read/write access to second proton 4-momentum
  FV& P2()  { return P[1]; }
  //! read only access to first proton 4-momentum
  FV const& P1() const { return P[0]; }
  //! read only access to second proton 4-momentum
  FV const& P2() const { return P[1]; }
  double const& get_rs()    const { return d_rs; }
  double get_s()            const { return std::pow(d_rs,2); }
  void set_rs(double const& rs) { d_rs=rs;}
  //! not used
  bool toggle_decay() { return d_decay = !d_decay; }
  //! get current phase space weight
  double const& get_wgt() { return d_wgt; }
  int set_parent(PS_2* parent) { if (parent != nullptr) { d_parent = parent; return 1; } return 0; }
  //! swap incoming parton and proton momenta
  void swap_initial_state();
  virtual double const& get_msq(int i) const = 0;
  virtual PS_2* get_child(int i) const = 0;
  virtual int set_child(int i, PS_2* child) = 0;
  virtual FV const& get_k(int i) const = 0;
  virtual int whattype() const = 0;
  virtual void print() const;

  virtual FV const& k1() const  { return nullvec; }
  virtual FV const& k2() const  { return nullvec; }
  virtual FV const& k3() const  { return nullvec; }
  /* virtual FV const& k1_lab() const  { return nullvec; } */
  /* virtual FV const& k2_lab() const  { return nullvec; } */
  /* virtual FV const& k3_lab() const  { return nullvec; } */
#ifdef WITH_T_SPIN  
  virtual FV const& s1() const  { return nullvec; }
  virtual FV const& s2() const  { return nullvec; }
  virtual FV const& s3() const  { return nullvec; }
  virtual FV const& s1_r() const  { return nullvec; }
  virtual FV const& s2_r() const  { return nullvec; }
  virtual FV const& s3_r() const  { return nullvec; }  
#endif
  
  static const FV nullvec;
};
////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
\typedef observable function type
 */
typedef double (*OBSFnc)(const PS_2*);

/*!
\typedef A distribution is specified by the histogram that stores the data and the observable that is historgammed.
 */
//typedef std::pair<HistArray*,OBSFnc> Distribution; 

/*!
\typedef A vector with pointers to distributions.
 */
typedef std::vector< std::shared_ptr<Distribution> > DistVec; 



////// class PS_2_1 ////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->1 scattering.
*/
class PS_2_1: public PS_2 {
 protected:
  FV k;
#ifdef WITH_T_SPIN   
  FV S;
#endif  
  double d_msq;
  PS_2* d_child;
  double d_x; // used in integrated dipoles: scaling factor for s
  
 public:
  explicit PS_2_1(const std::string& nm = "2->2");
  explicit PS_2_1(double const& msq, const std::string& nm = "2->2");
  ~PS_2_1();

  void   set_x(double const& x) { d_x = x;  } 
  double const& get_x()   const { return d_x;    }
  FV const& get_k(int i)  const { return k; }
  PS_2* get_child(int i)  const { return d_child; }
  
  double const& get_msq(int i)  const { return d_msq; }
  void set_msq(double const& msq) { d_msq = msq; }
  FV& k1()  { return k; }
  
  FV const& k1() const  { return k; }

  int set();
  
  int set_child(int i, PS_2* child) { if (child->set_parent(this)) { d_child = child; return 1;} return 0; }
  int whattype() const { return 21; }
  void print() const;
  void FillDistributions(
			 DistVec& dist,
			 H_TYPE id,
			 double const& wgt,
			 double const& mScale=1.0) const;
};
////// class PS_2_2 ////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->2 scattering.
*/
class PS_2_2: public PS_2 {
 private:

  //! outgoing particle 4-momenta
  FV k[2];
  /* //! outgoing particle 4-momenta in lab-frame */
  /* volatile FV k_lab[2];   */
#ifdef WITH_T_SPIN
  //! outgoing particle spin 4-vectors
  FV S[2];
  FV S_r[2];
#endif
  //! final state particle masses squared
  double d_msq[2];
  //! beta = P/E
  double d_beta[2];
  PS_2*  d_child[2];

  //! azimuthal scattering angle in the parton z.m.f.  (z-axis is the reference axis)
  double d_y;
  //! polar scattering angle in the parton z.m.f. (z-axis is the reference axis)
  double d_phi;
  //! t11 = 2 p1.k1
  double d_t11;
  //! t12 = 2 p1.k2
  double d_t12;
  //! initial state boost factor, used for integrated dipoles
  double d_x;
  //! beta*y
  double d_beta_y;

    
 public:

  explicit PS_2_2(const std::string& nm = "2->2");
  explicit PS_2_2(double const& m1sq, double const& m2sq, const std::string& nm = "2->2");
  ~PS_2_2();

  
  FV const& get_k(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return k[i]; }
  PS_2* get_child(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_child[i]; }
  double const& get_msq(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_msq[i]; }
  void set_msq(double const& m1sq,  double const& m2sq) { d_msq[0] = m1sq; d_msq[1] = m2sq; }
  double const& get_beta(int i) const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_beta[i]; }
  double const& get_beta()  const { return d_beta[0]; }
  double const& get_y()     const { return d_y;    }
  double const& get_phi()   const { return d_phi;  }
  double const& get_t11()   const { return d_t11;  }  
  double const& get_t12()   const { return d_t12;  }  
  double const& get_x()     const { return d_x;    }  
  double const& get_beta_y()const { return d_beta_y;    }  
  void   set_x(double const& x) { d_x = x;  }

  //! compute the phase space density for the current setting
  double cmp_wgt();
  
  // pointer to outgoing 4-momenta and spin
  FV& k1()  { return k[0]; }
  FV& k2()  { return k[1]; }
#ifdef WITH_T_SPIN    
  FV& s1()  { return S[0]; }
  FV& s2()  { return S[1]; }
  FV& s1_r()  { return S_r[0]; }
  FV& s2_r()  { return S_r[1]; }
#endif
  
  FV const& k1() const  { return k[0]; }
  FV const& k2() const  { return k[1]; }
  /* FV const& k1_lab() const  { return k_lab[0]; } */
  /* FV const& k2_lab() const  { return k_lab[1]; } */
#ifdef WITH_T_SPIN   
  FV const& s1() const  { return S[0]; }
  FV const& s2() const  { return S[1]; }
  FV const& s1_r() const  { return S_r[0]; }
  FV const& s2_r() const  { return S_r[1]; }
#endif
  
  //! boost initial state partons by a factor of d_x
  int boost_initial_state();
  //! boost initial state partons by a factor of x
  int boost_initial_state(double const& x);
  int boost_final_state();

  /*!
    New phase space setting: set 4-vectors according to c.m.e rs = sqrt(s_part) and scattering angles y, phi given in the parton z.m.f.
    \param rs \f$ \sqrt{s_{\rm part.}} \f$
    \param y azimuthal scattering angle in parton z.m.f. (z-axis is the reference axis)
    \param phi polar scattering angle in parton z.m.f. (z-axis is the reference axis)
  */
  int  set(double const& rs, 
	   double const& y, 
	   double const& phi=0.0);
  //! New phase space setting: compute c.m.e rs = sqrt(s_part) and scattering angles y, phi and all derived quantities from current set of 4-vectors. 
  int set(double const& x = 0.0);
  ////////////////////////////////////////

  int set_child(int i, PS_2* child) { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} if (child->set_parent(this)) { d_child[i] = child; return 1;} return 0; }
  int whattype() const { return 22; }
  void print() const;
  void FillDistributions(
			 DistVec& dist,
			 H_TYPE id,
			 double const& wgt,
			 double const& mScale=1.0) const;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_3 ////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->3 scattering.
*/
class PS_2_3: public PS_2 {
 private:
  //! outgoing particle 4-momenta
  FV k[3];
  /* //! outgoing particle 4-momenta in lab-frame */
  /* volatile FV k_lab[3]; */
#ifdef WITH_T_SPIN    
  //! outgoing particle spin 4-vectors
  FV S[3];
  FV S_r[3];
#endif
  
  //! final state particle masses squared
  double d_msq[3];
  //! beta = P/E
  double d_beta[3];
  PS_2*  d_child[3];

  double d_y_cm;
  double d_M12;
  double d_y_12;
  double d_phi_12;


 public:
  explicit PS_2_3(double const& m1sq, double const& m2sq, double const& m3sq, const std::string& nm = "2->3");
  explicit PS_2_3(const std::string& nm = "2->3");
  ~PS_2_3();

  // return the class members
  FV const& get_k(int i)  const { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} return k[i]; }
  PS_2* get_child(int i)  const { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} return d_child[i]; }
  double const& get_msq(int i)  const { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} return d_msq[i]; }
  void set_msq(double const& m1sq, double const& m2sq, double const& m3sq) { d_msq[0] = m1sq; d_msq[1] = m2sq; d_msq[2] = m3sq; }
  double const& get_beta(int i) const { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} return d_beta[i]; }
  double const& get_y_cm()  const { return d_y_cm;}
  double const& get_M12()   const { return d_M12;}
  double const& get_y_12()  const { return d_y_12;}
  double const& get_phi_12()const { return d_phi_12;}

  // pointer to outgoing 4-momenta and spin
  FV& k1()  { return k[0]; }
  FV& k2()  { return k[1]; }
  FV& k3()  { return k[2]; }
  // in the matrix elemets k3 is also called p3
  FV& p3()  { return k[2]; }
#ifdef WITH_T_SPIN    
  FV& s1()  { return S[0]; }
  FV& s2()  { return S[1]; }
  FV& s3()  { return S[2]; }
  FV& s1_r()  { return S_r[0]; }
  FV& s2_r()  { return S_r[1]; }
  FV& s3_r()  { return S_r[2]; }
#endif
  
  FV const& k1() const  { return k[0]; }
  FV const& k2() const  { return k[1]; }
  FV const& k3() const  { return k[2]; }
  FV const& p3() const  { return k[2]; }
  /* FV const& k1_lab() const  { return k_lab[0]; } */
  /* FV const& k2_lab() const  { return k_lab[1]; } */
  /* FV const& k3_lab() const  { return k_lab[2]; }   */
#ifdef WITH_T_SPIN    
  FV const& s1() const  { return S[0]; }
  FV const& s2() const  { return S[1]; }
  FV const& s3() const  { return S[2]; }
  FV const& s1_r() const  { return S_r[0]; }
  FV const& s2_r() const  { return S_r[1]; }
  FV const& s3_r() const  { return S_r[2]; }
#endif
  
  int boost_to_parent();

  /*!
    Compute new phase space setting: The 2->3 phase space is split into p1 p2 -> Q k3  and  Q -> k1 k2. y_cm is the angle of k3 in the p1 + p2 c.m.f., the mass of the intermediate state Q is M12 = sqrt((k1+k2)^2), the angles y_12 and phi_12 define the direction of k1' and k2' defined in the Q-restframe.
    \param rs \f$ \sqrt{s_{\rm part.}} \f$
    \param y_cm azimuthal scattering angle of k3 in p1+p2 z.m.f. (z-axis is the reference axis)
    \param phi_cm polar scattering angle of k3 in p1+p2 z.m.f. (z-axis is the reference axis)
    \param M12 invariant mass of the system k1+k2
    \param y_12  azimuthal scattering angle of k2 in k1+k2 z.m.f. (z-axis is the reference axis)
    \param phi_12  polar scattering angle of k2 in k1+k2 z.m.f. (z-axis is the reference axis)
  */
  int  set(double const& rs,
	   double const& y_cm,
	   double const& phi_cm,
	   double const& M12,
	   double const& y_12,
	   double const& phi_12);


  int set_child(int i, PS_2* child) { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} if (child->set_parent(this)) { d_child[i] = child; return 1;} return 0; }

  int whattype() const { return 23; }
  void print() const;
  void FillDistributions(
			 DistVec& dist,
			 H_TYPE id,
			 double const& wgt,
			 double const& mScale=1.0) const;
};
////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
// definition of observables used for this project, feel free to add more ///////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//! Computes invariant mass of two 4-vectors k1 and k2
double obs_M12(FV const& k1, FV const& k2);
double OBS_M12(const PS_2* ps);



//! Computes transverse momentum of a 4-vector k1 (transverse = x-y-plane), i.e. k_{T,1} = \sqrt{k1_1}^2+k1_{2}^2}
double obs_PT(FV const& k);
double OBS_PT1(const PS_2* ps);
double OBS_PT2(const PS_2* ps);


//! Computes transverse momentum of two particle system k1,k2, k_{T,12} = |k_{T,1}+k_{T,2}|
double obs_PT12(FV const& k1, FV const& k2);
double OBS_PT12(const PS_2* ps);


//! Computes rapidity of a 4-vector k1
double obs_Y(FV const& k);
double OBS_Y1(const PS_2* ps);
double OBS_Y2(const PS_2* ps);


//! Computes rapidity difference of two 4-vectors k1 and k2
double obs_DY(FV const& k1, FV const& k2);
double OBS_DY12(const PS_2* ps);



//! Computes spatial opening angle of two 4-vectors k1 and k2
double obs_PHI(FV const& k1, FV const& k2);

//! Computes opening angle of the spatial projections on the x-y plane
double obs_PHIT(FV const& k1, FV const& k2);

//! Computes the spatial triple product of three 4-vectors k1,k2,k3
double obs_TriProd(FV const& k1, FV const& k2, FV const& k3);


#ifdef WITH_T_SPIN  
double OBS_D12(const PS_2* ps);
double OBS_CP1(const PS_2* ps);
double OBS_CP2(const PS_2* ps);
double OBS_HEL12(const PS_2* ps);
double OBS_B1(const PS_2* ps);
double OBS_B2(const PS_2* ps);
#endif




#endif
