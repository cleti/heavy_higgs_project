
#ifndef SCALARINTEGRALS_H
#define SCALARINTEGRALS_H

#include "Global.h"
#include "Makros.h"
#include "Flags.h"
#include "PhaseSpace.h"

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

void Set_SI(const PS_2_2& ps, double MUR2, unsigned flags);
void Set_I1(double& MUR2);
void Set_I2(double& MUR2);
void Set_I3(double& MUR2);
void Set_I4(double& MUR2);
void Print_SI(int eps=0);



#endif
