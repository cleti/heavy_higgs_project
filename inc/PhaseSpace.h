
/*! \file
  \brief Phase space classes
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


////// class PS_2   ///////////////////////////////////////////////////////////////////////////////////////
/*!
  Abstract base class for 2->X phase space classes, holds among others the 4-vectors of the incoming partons.
  \sa Lorentz.h, HistArray.h
*/
class PS_2 {
 protected:
  //! square root of the current partonic c.m.e.
  double d_rs;
  //! name of the PS instance, useful to locate errors
  std::string d_name;
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
  explicit PS_2(const std::string& nm = "") : d_rs(0.0), d_name(nm), p{0.0,0.0}, P{0.0,0.0}, d_wgt(0.0), d_decay(false), d_parent(nullptr) {}
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
  double get_s()            const { return d_rs*d_rs; }
  void set_rs(double const& rs) { d_rs=rs;}
  void set_name(std::string const& name) { d_name = name; }
  //! not used
  bool toggle_decay() { return d_decay = !d_decay; }
  //! get current phase space weight
  double const& get_wgt() { return d_wgt; }
  int set_parent(PS_2* parent) { if (parent != nullptr) { d_parent = parent; return 1; } return 0; }
  //! swap incoming parton and proton momenta
  void swap();
  virtual double const& get_msq(int i) const = 0;
  virtual PS_2* get_child(int i) const = 0;
  virtual int set_child(int i, PS_2* child) = 0;
  virtual FV const& get_k(int i) const = 0;
  virtual int whattype() const = 0;
  virtual void print() const;
  virtual void FillDistributions(std::vector<HistArray*>& dist, int id, double const& wgt) const = 0;
};
////// class PS_2_1 ///////////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->1 scattering.
*/
class PS_2_1: public PS_2 {
 protected:
  FV k;
  FV S;
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
  void FillDistributions(std::vector<HistArray*>& dist, int id, double const& wgt) const;
};
////// class PS_2_2 ///////////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->2 scattering.
*/
class PS_2_2: public PS_2 {
 private:

  //! outgoing particle 4-momenta
  FV k[2];
  //! outgoing particle spin 4-vectors
  FV S[2];
  FV S_r[2];
  double d_msq[2];
  double d_beta[2];
  PS_2*  d_child[2];

  double d_y;
  double d_phi;
  double d_t11;
  double d_t12;
  double d_x; // used in integrated dipoles: scaling factor for s
  double d_beta_y;

    
 public:

  explicit PS_2_2(const std::string& nm = "2->2");
  explicit PS_2_2(double const& m1sq, double const& m2sq, const std::string& nm = "2->2");
  ~PS_2_2();

  
  // return the class members
  FV const& get_k(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return k[i]; }
  PS_2* get_child(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_child[i]; }
  double const& get_msq(int i)  const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_msq[i]; }
  void set_msq(double const& m1sq,  double const& m2sq) { d_msq[0] = m1sq; d_msq[1] = m2sq; }
  double const& get_beta(int i) const { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} return d_beta[i]; }
  double const& get_beta() const { return d_beta[0]; }
  double const& get_y()     const { return d_y;    }
  double const& get_phi()   const { return d_phi;  }
  double const& get_t11()   const { return d_t11;  }  // this is 2*p1.k1, not the invariant (p1-k1)^2 !
  double const& get_t12()   const { return d_t12;  }  // this is 2*p1.k2, not the invariant (p1-k2)^2 !
  double const& get_x()     const { return d_x;    }  // needed to evaluate integrated dipoles (
  double const& get_beta_y()const { return d_beta_y;    }  
  void   set_x(double const& x) { d_x = x;  } // needed to evaluate integrated dipoles (initial boost parameter)
  double cmp_wgt();
  
  // pointer to outgoing 4-momenta and spin
  FV& k1()  { return k[0]; }
  FV& k2()  { return k[1]; }
  FV& s1()  { return S[0]; }
  FV& s2()  { return S[1]; }
  FV& s1_r()  { return S_r[0]; }
  FV& s2_r()  { return S_r[1]; }

  FV const& k1() const  { return k[0]; }
  FV const& k2() const  { return k[1]; }
  FV const& s1() const  { return S[0]; }
  FV const& s2() const  { return S[1]; }
  FV const& s1_r() const  { return S_r[0]; }
  FV const& s2_r() const  { return S_r[1]; }
  
  // scale initial state vector by a factor of d_x
  int boost_initial_state();
  int boost_initial_state(double const& x);
  int boost_final_state();
  
  int  set(double const& rs, 
	   double const& y, 
	   double const& phi=0.0);
  void set();
  ////////////////////////////////////////

  int set_child(int i, PS_2* child) { if (!RANGE(i,2)) {ERROR("i is out of range (2)");} if (child->set_parent(this)) { d_child[i] = child; return 1;} return 0; }
  int whattype() const { return 22; }
  void print() const;
  void FillDistributions(std::vector<HistArray*>& dist, int id, double const& wgt) const;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////
////// class PS_2_3 ///////////////////////////////////////////////////////////////////////////////////////
/*!
  Phase space for 2->3 scattering.
*/
class PS_2_3: public PS_2 {
 private:
  //! outgoing particle 4-momenta
  FV k[3];
  //! outgoing particle spin 4-vectors
  FV S[3];
  FV S_r[3];
  
  double d_msq[3];
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
  FV& s1()  { return S[0]; }
  FV& s2()  { return S[1]; }
  FV& s1_r()  { return S_r[0]; }
  FV& s2_r()  { return S_r[1]; }
  
  FV const& k1() const  { return k[0]; }
  FV const& k2() const  { return k[1]; }
  FV const& k3() const  { return k[2]; }
  FV const& p3() const  { return k[2]; }
  FV const& s1() const  { return S[0]; }
  FV const& s2() const  { return S[1]; }
  FV const& s1_r() const  { return S_r[0]; }
  FV const& s2_r() const  { return S_r[1]; }
  
  int boost_to_parent();
  // 2->3 phase space is split into p1 p2 -> Q k3  and  Q -> k1 k2
  // y_cm is the angle of k3 in the p1 + p2 c.m.f., the mass of the intermediate state M12 = sqrt(Q^2)
  // the angles y_12 and phi_12 define the direction of k1' and k2' defined in the Q-restframe
  // default masses: k1^2 = 1, k2^2 = 1, k3^2 = 0 corresponding to tt+g production
  int  set(double const& rs,
	   double const& y_cm,
	   double const& phi_cm,
	   double const& M12,
	   double const& y_12,
	   double const& phi_12);


  int set_child(int i, PS_2* child) { if (!RANGE(i,3)) {ERROR("i is out of range (3)");} if (child->set_parent(this)) { d_child[i] = child; return 1;} return 0; }

  int whattype() const { return 23; }
  void print() const;
  void FillDistributions(std::vector<HistArray*>& dist, int id, double const& wgt) const;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////








#endif
