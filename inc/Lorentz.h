
/*! \file
  \brief This file provides the 4-vector and Lorentz transformation classes. 
*/

#ifndef LORENTZ_H
#define LORENTZ_H

#include <initializer_list>

#include "Makros.h"
#include <vector>
#include <iostream>

////// class FV     ///////////////////////////////////////////////////////////////////////////////////////
#define RANGE(i,J) (i>=0 && i<J)
/*!
  \brief Lorentz 4-vector.
*/
class FV {
 protected:
  //! the 4-vector components
  std::vector<double> v;//double v[4];

 public:
  FV(double const& a=0.0)                      : v(4,a)              { }
  FV(FV const& rhs)                            : v(rhs.v)            { }
  FV(FV&& rhs) noexcept                        : v(std::move(rhs.v)) { }
  FV(std::initializer_list<double> const& rhs) : v(rhs)              { }
  FV(std::vector<double> const&  vrhs)         : v(vrhs)             { }
  ~FV() { /* delete v; */ }


  /* FV(double const& a=0.0)                   : v{a,a,a,a} {    } */
  /* FV(FV const& rhs)                         : v(rhs.v)   {  }  { *this = rhs} */
  /* FV(FV&& rhs) noexcept /\*                 : FV()  *\/  { *this=std::move(rhs); } */
  /* FV(std::initializer_list<double> rhs) /\* : FV() *\/   { *this=rhs; } */
  /* ~FV() { /\* delete v; *\/ } */

  //! Read/write access to specific component. NO RANGE CHECK!
  double& operator[](int const& i)             {  return v[i]; }
  //! Read access to specific component. NO RANGE CHECK!
  double const& operator[](int const& i) const {  return v[i]; }
  FV& operator=(std::initializer_list<double> L);
  FV& operator=(FV const& other);
  FV& operator=(FV&& other) noexcept;
  FV& operator+=(FV const& other);
  FV& operator-=(FV const& other);
  FV& operator*=(double const& a);
  FV& operator/=(double const& a);

  //! Swap vectors, using std::swap(std::vector) internally.
  void swap(FV& other);
  //! Normalize spatial components to one. Returns a copy.
  FV norm3() const;
};

FV operator*(FV const& v ,double const& a);
FV operator*(double const& a, FV const& v);
FV&& operator*(FV&& v ,double const& a);
FV&& operator*(double const& a, FV&& v);

FV operator/(FV const& v ,double const& a);
FV&& operator/(FV&& v ,double const& a);

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
/*!
  \brief Lorentz transformation.
*/
class LT {
 protected:
  //! 4 x 4 components of the matrix.
  double M[4][4];
  //! The metric tensor G = [+1,-1,-1,-1].
  static const double G[4];

 public:
  LT();
  ~LT();

  /* int set_col(int i, double* vals); */
  /* int set_row(int i, double* vals); */
  //! Invert the Lorentz transformation. This does not compute the inverse of general 4 x 4 matrices!
  void invert();
  //! Transpose matrix.
  void transpose();
  /*!
    Apply transformation to covariant components of the 4-vector.
    \param v 4-vector to be transformed.
  */
  void apply(FV& v);
  /*!
    Apply transformation to contravariant components of the 4-vector src and store result in trg.
    \param src source 4-vector.
    \param trg target 4-vector.
  */  
  void apply_cpy(FV const &src, FV &trg);
    
  /*!
    Apply transformation to contravariant components of the 4-vector, i.e. multiplication with metric tensor G.
    \param v 4-vector to be transformed.
  */
  void apply_G(FV& v);
  /*!
    Apply transformation to contravariant components of the 4-vector src and store result in trg.
    \param src source 4-vector.
    \param trg target 4-vector.
  */  
  void apply_G_cpy(FV const &src, FV &trg);
    
  int set_FF(FV const& p1, FV const& p2);
  int set_II(FV const& K , FV const& Kb);
  /*!
    Set boost that relates the frame in which P is given to its restframe.
    \param P defines the boost
    \param INV invert boost, same as set_boost_inv()
    \returns 1 if succesful, 0 if the mass of P is too small so that the restframe does not exist.
  */
  int set_boost(FV const& P,bool INV=false);
  /*!
    Set inverse boost that relates the P restframe to the frame in which P is given.
    \param defines the boost
    \returns 1 if succesful, 0 if the mass of P is too small so that the restframe does not exist.
  */
  int set_boost_inv(FV const& P);
  
  int set_boost_z(double const& x,bool INV=false);
  
  int set_wigner(FV const& P1, FV const& P2);
  void print();
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////




#endif
