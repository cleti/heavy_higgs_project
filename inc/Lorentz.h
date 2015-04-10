

#ifndef LORENTZ_H
#define LORENTZ_H

#include <initializer_list>

#include "Makros.h"


////// class FV     ///////////////////////////////////////////////////////////////////////////////////////
#define RANGE(i,J) (i>=0 && i<J)

class FV {
 protected:
  /* double* v; */
  double v[4];

 public:
  //FV(double const& a=0.0) : v(new double[4])    { if (a!=0.0) {v[0]=a;v[1]=a;v[2]=a;v[3]=a;} }std::cout << std::endl << " std.  constr. " << std::endl;
  FV(double const& a=0.0) : v{a,a,a,a}       {    }
  FV(FV const& rhs) : v{rhs[0],rhs[1],rhs[2],rhs[3]} {  }
  FV(FV&& rhs) noexcept /* : FV()  */                 { *this=std::move(rhs); }
  FV(std::initializer_list<double> rhs) /* : FV() */  { *this=rhs; }
  ~FV() { /* delete v; */ }

  // BEWARE! no range check here!!!
  //if (!RANGE(i,4)) {ERROR("out of range");}
  double& operator[](int const& i)             {  return v[i]; }
  double const& operator[](int const& i) const {  return v[i]; }
  FV& operator=(std::initializer_list<double> L);
  FV& operator=(FV const& other);
  FV& operator=(FV&& other) noexcept;
  FV& operator+=(FV const& other);
  FV& operator-=(FV const& other);
  FV& operator*=(double const& a);
  FV& operator/=(double const& a);
};

FV operator*(FV const& v ,double const& a);
FV operator*(double const& a, FV const& v);

FV operator+(FV const& v1,FV const& v2);
FV&& operator+(FV const& v1,FV&& v2) noexcept;
FV&& operator+(FV&& v1,FV const& v2) noexcept;
FV&& operator+(FV&& v1,FV&& v2) noexcept;

FV operator-(FV const& v1,FV const& v2);
FV&& operator-(FV const& v1,FV&& v2) noexcept;
FV&& operator-(FV&& v1,FV const& v2) noexcept;
FV&& operator-(FV&& v1,FV&& v2) noexcept;

FV operator-(FV const& v1);
FV&& operator-(FV&& v1) noexcept;





////// class LT     ///////////////////////////////////////////////////////////////////////////////////////
class LT {
 protected:
  double M[4][4];
  static const double G[4];

 public:
  LT();
  ~LT();

  /* int set_col(int i, double* vals); */
  /* int set_row(int i, double* vals); */
  void invert();
  void transpose();
  void apply(FV& v);
  void apply_G(FV& v);
  int set_FF(FV const& p1, FV const& p2);
  int set_II(FV const& K , FV const& Kb);
  // boost vectors defined in the current frame (where also the argument P is defined) to P-restframe
  int set_boost(FV const& P,bool INV=false);
  // boost vectors from P-restframe to the current frame (where also the argument P is defined)
  int set_boost_inv(FV const& P);

  int set_wigner(FV const& P1, FV const& P2);
  void print();
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////




#endif
