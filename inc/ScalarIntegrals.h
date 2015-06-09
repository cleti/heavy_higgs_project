

/*! \file
  \brief This file provides wrapper functions for the scalar integrals provided by the QCDloop  library arXiv:0712.1851 [hep-ph]. A number of global variables is defined which is used for the evaluation of the one-loop amplitudes relevant for this project. 
*/


#ifndef SCALARINTEGRALS_H
#define SCALARINTEGRALS_H

#include "Global.h"
#include "Makros.h"
#include "Flags.h"

#ifdef WITH_NON_FACT_DIAGRAMS
#include "../LoopTools-2.12/build/clooptools.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif
  extern c_double qli1_(const double*,const double*,const int*);
  extern c_double qli2_(const double*,const double*,const double*,const double*,const int*);
  extern c_double qli3_(const double*,const double*,const double*,const double*,const double*,const double*,const double*,const int*);
  extern c_double qli4_(const double*,const double*,const double*,const double*,const double*,const double*,const double*,const double*,const double*,const double*,const double*,const int*);
#ifdef __cplusplus
}
#endif


extern c_double I1_MT2_MU2_0;
/* extern c_double I2_0_MT2_MT2_MU2_0; */
extern c_double I2_MT2_0_MT2_MU2_0;
extern c_double I2_S12_0_0_MU2_0;
extern c_double I2_S12_MT2_MT2_MU2_0; 
extern c_double I2_T11_0_MT2_MU2_0;
extern c_double I2_T12_0_MT2_MU2_0;
extern c_double I3_0_0_S12_0_0_0_MU2_0;
extern c_double I3_0_T11_MT2_0_0_MT2_MU2_0;
extern c_double I3_0_T12_MT2_0_0_MT2_MU2_0;
extern c_double I3_MT2_0_T11_0_MT2_MT2_MU2_0;
extern c_double I3_MT2_0_T12_0_MT2_MT2_MU2_0;
extern c_double I3_MT2_S12_MT2_0_MT2_MT2_MU2_0;
extern c_double I3_S12_0_0_MT2_MT2_MT2_MU2_0;
extern c_double I3_S12_MT2_MT2_0_0_MT2_MU2_0;
extern c_double I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_0;
extern c_double I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_0;
extern c_double I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_0;
extern c_double I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_0;
extern c_double I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0;
extern c_double I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0;

extern c_double I1_MT2_MU2_1;
/* extern c_double I2_0_MT2_MT2_MU2_1; */
extern c_double I2_MT2_0_MT2_MU2_1;
extern c_double I2_S12_0_0_MU2_1;
/* extern c_double I2_S12_MT2_MT2_MU2_1; */
extern c_double I2_T11_0_MT2_MU2_1;
extern c_double I2_T12_0_MT2_MU2_1;
extern c_double I3_0_0_S12_0_0_0_MU2_1;
extern c_double I3_0_T11_MT2_0_0_MT2_MU2_1;
extern c_double I3_0_T12_MT2_0_0_MT2_MU2_1;
/* extern c_double I3_MT2_0_T11_0_MT2_MT2_MU2_1; */
/* extern c_double I3_MT2_0_T12_0_MT2_MT2_MU2_1; */
extern c_double I3_MT2_S12_MT2_0_MT2_MT2_MU2_1;
/* extern c_double I3_S12_0_0_MT2_MT2_MT2_MU2_1; */
/* extern c_double I3_S12_MT2_MT2_0_0_MT2_MU2_1; */
extern c_double I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_1;
extern c_double I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_1;
extern c_double I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_1;
extern c_double I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_1;
extern c_double I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1;
extern c_double I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1;


#ifdef WITH_NON_FACT_DIAGRAMS
// scalar integrals from the non-fact. loops
extern c_double I1_MH2_MU2_0;
extern c_double I1_MH2_MU2_1;

extern c_double I2_S12_0_MH2_MU2_0;
extern c_double I2_MT2_MT2_MH2_MU2_0;
extern c_double I2_0_0_MH2_MU2_0;
extern c_double I2_S12_0_MH2_MU2_1;
extern c_double I2_MT2_MT2_MH2_MU2_1;
extern c_double I2_0_0_MH2_MU2_1;

extern c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_0;
extern c_double I3_T11_MT2_0_0_MT2_MH2_MU2_0;
extern c_double I3_T12_MT2_0_0_MT2_MH2_MU2_0;
extern c_double I3_0_0_S12_0_0_MH2_MU2_0;
extern c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_1;
extern c_double I3_T11_MT2_0_0_MT2_MH2_MU2_1;
extern c_double I3_T12_MT2_0_0_MT2_MH2_MU2_1;
extern c_double I3_0_0_S12_0_0_MH2_MU2_1;

extern c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_2;
extern c_double I3_T11_MT2_0_0_MT2_MH2_MU2_2;
extern c_double I3_T12_MT2_0_0_MT2_MH2_MU2_2;
extern c_double I3_0_0_S12_0_0_MH2_MU2_2;
#endif

/*!
  \brief Set global variables according to the given phase space point.
  \param mt2  top-quark mass squared
  \param rs   sqrt(s)
  \param t11  2*p1.k1
  \param t12  2*p1.k2
  \param MUR2 renormalization scale squared
  \param flags specify which groups of integrals get evaluated
  \sa Flags.h
*/
void Set_SI(
	    const double& mt2,
	    const double& rs,
	    const double& t11,
	    const double& t12,
	    const double& MUR2,
	    unsigned flags);
/*!
  \brief Evaluate the 1-point functions.
*/
void Set_I1(const double& MUR2);
/*! 
  \brief Evaluate the 2-point functions.
*/
void Set_I2(const double& MUR2);
/*! 
  \brief Evaluate the 3-point functions.
*/
void Set_I3(const double& MUR2);
/*!
  \brief Evaluate the 4-point functions.
*/
void Set_I4(const double& MUR2);
void Print_SI(int eps=0);



#endif
