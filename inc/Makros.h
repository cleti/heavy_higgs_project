
/*! \file
  \brief Makros for printing error, warning messages, etc.
 */


#ifndef MAKROS_H
#define MAKROS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <thread>
#include <boost/utility.hpp>




//! Print warning to stream OST.
#define SWARNING(OST,MSG)			\
  {						\
    OST << std::endl;				\
    OST << " Warning in ";			\
    OST << __FILE__ << ", ";			\
    OST << __FUNCTION__  << ", line ";		\
    OST << __LINE__ << ": " <<  MSG << ". ";	\
    OST << std::endl;				\
  }						\

#define WARNING(MSG)	SWARNING(std::cout,MSG)


//! Print error message to stream OST without exit.
#define SERROR_NOEXIT(OST,MSG)		       	\
  {							\
    OST << std::endl;					\
    std::cout << " ERROR in ";				\
    std::cout << __FILE__ << ", ";			\
    std::cout << __FUNCTION__  << ", line ";		\
    std::cout << __LINE__ << ": " <<  MSG << "! ";	\
    std::cout << std::endl;				\
  }							\

//! Print error message to stream OST and exit.
#define SERROR(OST,MSG)		       	\
  {						\
    ERROR_NOEXIT(MSG);				\
    std::exit(1);				\
  }						\

//! Print error message to standard output without exit.
#define ERROR_NOEXIT(MSG)   SERROR_NOEXIT(std::cout,MSG);

//! Print error message to standard output and exit.
#define ERROR(MSG)  SERROR(std::cout,MSG);



//! check float for inf/nan, print ps content
#define CHECKNAN_PS(VAR,PS)			\
  if ( std::isnan(VAR) || std::isinf(VAR) )	\
    {						\
      ERROR_NOEXIT(#VAR << " is inf/nan");	\
      PS.print();				\
      std::exit(1);				\
    }						\

//! check float for inf/nan
#define CHECKNAN_NOEXIT(VAR)			\
  if ( std::isnan(VAR) || std::isinf(VAR) )	\
    {						\
      ERROR_NOEXIT(#VAR << " is inf/nan");	\
    }						\

//! check float for inf/nan and exit
#define CHECKNAN(VAR)				\
  if ( std::isnan(VAR) || std::isinf(VAR) )	\
    {						\
      ERROR(#VAR << " is inf/nan");		\
    }						\


#define unlikely(expr) __builtin_expect(!!(expr), 0)
#define likely(expr) __builtin_expect(!!(expr), 1)

//! delete a pointer and reset
#define DELETE_PTR(ptr)		      	\
  {						\
    delete ptr;					\
    ptr = nullptr;				\
  }						\




#define SLEEP(S)						\
  {								\
    std::this_thread::sleep_for( std::chrono::seconds(S) );	\
  }								\


#define EXIT(S)			       	\
  {						\
    std::cout << std::endl;			\
    std::cout << " exit called from ";		\
    std::cout << __FILE__ << ", ";		\
    std::cout << __FUNCTION__ << ", line ";	\
    std::cout << __LINE__ << std::endl;		\
    exit(S);					\
  }						\


//! Print variable name and value to ouput stream.
#define SPRINT(OST,VAR)		       	\
  {						\
    OST << std::endl;				\
    OST << " Value of ";			\
    OST << std::setw(20) << #VAR << " = ";	\
    OST << std::setw(20) << VAR;		\
    OST << std::endl;				\
  }						\

//! Print variable name and value to standard output.
#define PRINT(VAR)   SPRINT(std::cout,VAR)



#define RE(X) X.real()
#define IM(X) X.imag()
#define VF(X) (X)


//! Lorentz product of two 4 vectors K1, K2.
#define sp(K1,K2) (K1[0]*K2[0]-(K1[1]*K2[1]+K1[2]*K2[2]+K1[3]*K2[3]))

//! Scalar product of the spatial components of 4-vectors K1, K2.
#define sp3(K1,K2) (K1[1]*K2[1]+K1[2]*K2[2]+K1[3]*K2[3])

//! Mass squared of a 4-vector K.
#define MSQ(K) sp(K,K)

//! Length of a 4-vector K.
#define LEN(K) std::sqrt(K[1]*K[1]+K[2]*K[2]+K[3]*K[3])

//! Print 4-vector P.
#define PRINT_4VEC(P)							\
  {									\
    std::cout << std::endl << #P << " =  [ ";				\
      std::cout << std::setw(8) << P[0] << ", ";			\
      std::cout << std::setw(8) << P[1] << ", ";			\
      std::cout << std::setw(8) << P[2] << ", ";			\
      std::cout << std::setw(8) << P[3] << " ]";			\
      std::cout << std::endl << #P << "^2 = " << std::setw(8) << MSQ(P) << std::endl; \
  }									\

//! We use the GSL dilogarithm function in this program.
#define dilog(x) gsl_sf_dilog(x)



#endif
