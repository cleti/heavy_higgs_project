

#include "../inc/ScalarIntegrals.h"

// these are exported
c_double I1_MT2_MU2_0 = 0.0;
// c_double I2_0_MT2_MT2_MU2_0 = 0.0;
c_double I2_MT2_0_MT2_MU2_0 = 0.0;
c_double I2_S12_0_0_MU2_0 = 0.0;
c_double I2_S12_MT2_MT2_MU2_0 = 0.0;
c_double I2_T11_0_MT2_MU2_0 = 0.0;
c_double I2_T12_0_MT2_MU2_0 = 0.0;

c_double I3_0_0_S12_0_0_0_MU2_0 = 0.0;
c_double I3_0_T11_MT2_0_0_MT2_MU2_0 = 0.0;
c_double I3_0_T12_MT2_0_0_MT2_MU2_0 = 0.0;
c_double I3_MT2_0_T11_0_MT2_MT2_MU2_0 = 0.0;
c_double I3_MT2_0_T12_0_MT2_MT2_MU2_0 = 0.0;
c_double I3_MT2_S12_MT2_0_MT2_MT2_MU2_0 = 0.0;
c_double I3_S12_0_0_MT2_MT2_MT2_MU2_0 = 0.0;
c_double I3_S12_MT2_MT2_0_0_MT2_MU2_0 = 0.0;
c_double I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_0 = 0.0;
c_double I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_0 = 0.0;
c_double I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_0 = 0.0;
c_double I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_0 = 0.0;
c_double I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0 = 0.0;
c_double I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0 = 0.0;

c_double I1_MT2_MU2_1 = 0.0;
// c_double I2_0_MT2_MT2_MU2_1 = 0.0;
c_double I2_MT2_0_MT2_MU2_1 = 0.0;
c_double I2_S12_0_0_MU2_1 = 0.0;
// c_double I2_S12_MT2_MT2_MU2_1 = 0.0;
c_double I2_T11_0_MT2_MU2_1 = 0.0;
c_double I2_T12_0_MT2_MU2_1 = 0.0;
c_double I3_0_0_S12_0_0_0_MU2_1 = 0.0;
c_double I3_0_T11_MT2_0_0_MT2_MU2_1 = 0.0;
c_double I3_0_T12_MT2_0_0_MT2_MU2_1 = 0.0;
// c_double I3_MT2_0_T11_0_MT2_MT2_MU2_1 = 0.0;//
// c_double I3_MT2_0_T12_0_MT2_MT2_MU2_1 = 0.0;//
c_double I3_MT2_S12_MT2_0_MT2_MT2_MU2_1 = 0.0;
// c_double I3_S12_0_0_MT2_MT2_MT2_MU2_1 = 0.0;//
// c_double I3_S12_MT2_MT2_0_0_MT2_MU2_1 = 0.0;//
c_double I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_1 = 0.0;
c_double I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_1 = 0.0;
c_double I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_1 = 0.0;
c_double I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_1 = 0.0;
c_double I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1 = 0.0;
c_double I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1 = 0.0;


#ifdef WITH_NON_FACT_DIAGRAMS
// scalar integrals from the non-fact. loops
c_double I1_MH2_MU2_0 = 0.0;
c_double I1_MH2_MU2_1 = 0.0;

c_double I2_S12_0_MH2_MU2_0 = 0.0;
c_double I2_MT2_MT2_MH2_MU2_0 = 0.0;
c_double I2_0_0_MH2_MU2_0 = 0.0;

c_double I2_S12_0_MH2_MU2_1 = 0.0;
c_double I2_MT2_MT2_MH2_MU2_1 = 0.0;
c_double I2_0_0_MH2_MU2_1 = 0.0;


c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_0 = 0.0;
c_double I3_T11_MT2_0_0_MT2_MH2_MU2_0 = 0.0;
c_double I3_T12_MT2_0_0_MT2_MH2_MU2_0 = 0.0;
c_double I3_0_0_S12_0_0_MH2_MU2_0 = 0.0;

c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_1 = 0.0;
c_double I3_T11_MT2_0_0_MT2_MH2_MU2_1 = 0.0;
c_double I3_T12_MT2_0_0_MT2_MH2_MU2_1 = 0.0;
c_double I3_0_0_S12_0_0_MH2_MU2_1 = 0.0;

c_double I3_MT2_MT2_S12_0_MT2_MH2_MU2_2 = 0.0;
c_double I3_T11_MT2_0_0_MT2_MH2_MU2_2 = 0.0;
c_double I3_T12_MT2_0_0_MT2_MH2_MU2_2 = 0.0;
c_double I3_0_0_S12_0_0_MH2_MU2_2 = 0.0;
#endif


namespace SI_INV
{
  double S12 = 0.0;
  double T11 = 0.0;
  double T12 = 0.0;
  double zero = 0.0;
  c_double MH2 = 0.0;
}

void Set_SI(const PS_2_2& ps, double MUR2, unsigned flags)
{
  using namespace SI_INV;
  using RunParameters::mt2;
 
#ifdef WITH_NON_FACT_DIAGRAMS
  setmudim(MUR2);
  MH2 = pow(c_double(Bosons::M_1,Bosons::G_1),2);
  PRINT(MH2);
#endif
  
if (IZ_EVAL_SI_2P(flags))
  {

    // QCDloop functions take these invariants as arguments
    S12  = pow(ps.get_rs(),2);
    T11  = mt2 - ps.get_t11();
    T12  = mt2 - ps.get_t12();

    Set_I1(MUR2);
    Set_I2(MUR2);
    
    if (IZ_EVAL_SI_3P(flags))
      {
	Set_I3(MUR2);
	if (IZ_EVAL_SI_4P(flags))
	  {
	    Set_I4(MUR2);
	  }
      }
  }
}


void Set_I1(double& MUR2)
{
  using namespace SI_INV;
  using RunParameters::mt2;
    
  int I = 0;
  I1_MT2_MU2_0 = qli1_(&mt2,&MUR2,&I);
  I = -1;
  I1_MT2_MU2_1 = qli1_(&mt2,&MUR2,&I);

#ifdef WITH_NON_FACT_DIAGRAMS
  setlambda(0);
  I1_MH2_MU2_0 = A0C(MH2);
  setlambda(-1);
  I1_MH2_MU2_1 = A0C(MH2);
#endif

}

void Set_I2(double& MUR2)
{
  using namespace SI_INV;
  using RunParameters::mt2;
    
  int I = 0;
  // I2_0_MT2_MT2_MU2_0   = qli2_(&zero,&mt2,&mt2,&MUR2,&I);
  I2_MT2_0_MT2_MU2_0   = qli2_(&mt2,&zero,&mt2,&MUR2,&I);
  I2_S12_0_0_MU2_0     = qli2_(&S12,&zero,&zero,&MUR2,&I);
  I2_S12_MT2_MT2_MU2_0 = qli2_(&S12,&mt2,&mt2,&MUR2,&I);
  I2_T11_0_MT2_MU2_0   = qli2_(&T11,&zero,&mt2,&MUR2,&I);
  I2_T12_0_MT2_MU2_0   = qli2_(&T12,&zero,&mt2,&MUR2,&I);
 
  I = -1;
  // I2_0_MT2_MT2_MU2_1   = qli2_(&zero,&mt2,&mt2,&MUR2,&I);
  I2_MT2_0_MT2_MU2_1   = qli2_(&mt2,&zero,&mt2,&MUR2,&I);
  I2_S12_0_0_MU2_1     = qli2_(&S12,&zero,&zero,&MUR2,&I);
  // I2_S12_MT2_MT2_MU2_1 = qli2_(&S12,&mt2,&mt2,&MUR2,&I);
  I2_T11_0_MT2_MU2_1   = qli2_(&T11,&zero,&mt2,&MUR2,&I);
  I2_T12_0_MT2_MU2_1   = qli2_(&T12,&zero,&mt2,&MUR2,&I);


#ifdef WITH_NON_FACT_DIAGRAMS
  setlambda(0);
  I2_S12_0_MH2_MU2_0   = B0C(S12,zero,MH2);
  I2_MT2_MT2_MH2_MU2_0 = B0C(mt2,mt2,MH2);
  I2_0_0_MH2_MU2_0     = B0C(zero,zero,MH2);
  setlambda(-1);
  I2_S12_0_MH2_MU2_1   = B0C(S12,zero,MH2);
  I2_MT2_MT2_MH2_MU2_1 = B0C(mt2,mt2,MH2);
  I2_0_0_MH2_MU2_1     = B0C(zero,zero,MH2);
#endif
}

 void Set_I3(double& MUR2)
 {
  using namespace SI_INV;
  using RunParameters::mt2;
    
  int I = 0;
  I3_0_0_S12_0_0_0_MU2_0         = qli3_(&zero,&zero,&S12,&zero,&zero,&zero,&MUR2,&I);
  I3_0_T11_MT2_0_0_MT2_MU2_0     = qli3_(&zero,&T11,&mt2,&zero,&zero,&mt2,&MUR2,&I);
  I3_0_T12_MT2_0_0_MT2_MU2_0     = qli3_(&zero,&T12,&mt2,&zero,&zero,&mt2,&MUR2,&I);
  I3_MT2_0_T11_0_MT2_MT2_MU2_0   = qli3_(&mt2,&zero,&T11,&zero,&mt2,&mt2,&MUR2,&I);
  I3_MT2_0_T12_0_MT2_MT2_MU2_0   = qli3_(&mt2,&zero,&T12,&zero,&mt2,&mt2,&MUR2,&I);
  I3_MT2_S12_MT2_0_MT2_MT2_MU2_0 = qli3_(&mt2,&S12,&mt2,&zero,&mt2,&mt2,&MUR2,&I);
  I3_S12_0_0_MT2_MT2_MT2_MU2_0   = qli3_(&S12,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I);
  I3_S12_MT2_MT2_0_0_MT2_MU2_0   = qli3_(&S12,&mt2,&mt2,&zero,&zero,&mt2,&MUR2,&I);

  I = -1;
  I3_0_0_S12_0_0_0_MU2_1         = qli3_(&zero,&zero,&S12,&zero,&zero,&zero,&MUR2,&I);
  I3_0_T11_MT2_0_0_MT2_MU2_1     = qli3_(&zero,&T11,&mt2,&zero,&zero,&mt2,&MUR2,&I);
  I3_0_T12_MT2_0_0_MT2_MU2_1     = qli3_(&zero,&T12,&mt2,&zero,&zero,&mt2,&MUR2,&I);
  // I3_MT2_0_T11_0_MT2_MT2_MU2_1   = qli3_(&mt2,&zero,&T11,&zero,&mt2,&mt2,&MUR2,&I);//
  // I3_MT2_0_T12_0_MT2_MT2_MU2_1   = qli3_(&mt2,&zero,&T12,&zero,&mt2,&mt2,&MUR2,&I);//
  I3_MT2_S12_MT2_0_MT2_MT2_MU2_1 = qli3_(&mt2,&S12,&mt2,&zero,&mt2,&mt2,&MUR2,&I);
  // I3_S12_0_0_MT2_MT2_MT2_MU2_1   = qli3_(&S12,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I);//
  // I3_S12_MT2_MT2_0_0_MT2_MU2_1   = qli3_(&S12,&mt2,&mt2,&zero,&zero,&mt2,&MUR2,&I);//

  
#ifdef WITH_NON_FACT_DIAGRAMS
  setlambda(0);
  I3_MT2_MT2_S12_0_MT2_MH2_MU2_0 = C0C(mt2,mt2,S12,zero,mt2,MH2);
  I3_T11_MT2_0_0_MT2_MH2_MU2_0   = C0C(T11,mt2,zero,zero,mt2,MH2);
  I3_T12_MT2_0_0_MT2_MH2_MU2_0   = C0C(T12,mt2,zero,zero,mt2,MH2);
  I3_0_0_S12_0_0_MH2_MU2_0       = C0C(zero,zero,S12,zero,zero,MH2);  
  setlambda(-1);
  I3_MT2_MT2_S12_0_MT2_MH2_MU2_1 = C0C(mt2,mt2,S12,zero,mt2,MH2);
  I3_T11_MT2_0_0_MT2_MH2_MU2_1   = C0C(T11,mt2,zero,zero,mt2,MH2);
  I3_T12_MT2_0_0_MT2_MH2_MU2_1   = C0C(T12,mt2,zero,zero,mt2,MH2);
  I3_0_0_S12_0_0_MH2_MU2_1       = C0C(zero,zero,S12,zero,zero,MH2);
  setlambda(-2);
  I3_MT2_MT2_S12_0_MT2_MH2_MU2_2 = C0C(mt2,mt2,S12,zero,mt2,MH2);
  I3_T11_MT2_0_0_MT2_MH2_MU2_2   = C0C(T11,mt2,zero,zero,mt2,MH2);
  I3_T12_MT2_0_0_MT2_MH2_MU2_2   = C0C(T12,mt2,zero,zero,mt2,MH2);
  I3_0_0_S12_0_0_MH2_MU2_2       = C0C(zero,zero,S12,zero,zero,MH2);
#endif
}

 void Set_I4(double& MUR2)
 {
  using namespace SI_INV;
  using RunParameters::mt2;
    
  int I = 0;
  I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0 = qli4_(&mt2,&zero,&zero,&mt2,&T11,&S12,&zero,&mt2,&mt2,&mt2,&MUR2,&I);
  I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0 = qli4_(&mt2,&zero,&zero,&mt2,&T12,&S12,&zero,&mt2,&mt2,&mt2,&MUR2,&I);

  I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_0     = qli4_(&zero,&zero,&mt2,&mt2,&S12,&T11,&zero,&zero,&zero,&mt2,&MUR2,&I);
  I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_0     = qli4_(&zero,&zero,&mt2,&mt2,&S12,&T12,&zero,&zero,&zero,&mt2,&MUR2,&I);

  // these two are equal
  I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_0   = qli4_(&zero,&mt2,&zero,&mt2,&T11,&T12,&zero,&zero,&mt2,&mt2,&MUR2,&I);
  I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_0   = qli4_(&zero,&mt2,&zero,&mt2,&T12,&T11,&zero,&zero,&mt2,&mt2,&MUR2,&I);


  I = -1;
  I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1 = qli4_(&mt2,&zero,&zero,&mt2,&T11,&S12,&zero,&mt2,&mt2,&mt2,&MUR2,&I);
  I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1 = qli4_(&mt2,&zero,&zero,&mt2,&T12,&S12,&zero,&mt2,&mt2,&mt2,&MUR2,&I);

  I4_0_0_MT2_MT2_S12_T11_0_0_0_MT2_MU2_1     = qli4_(&zero,&zero,&mt2,&mt2,&S12,&T11,&zero,&zero,&zero,&mt2,&MUR2,&I);
  I4_0_0_MT2_MT2_S12_T12_0_0_0_MT2_MU2_1     = qli4_(&zero,&zero,&mt2,&mt2,&S12,&T12,&zero,&zero,&zero,&mt2,&MUR2,&I);

  I4_0_T11_0_T12_MT2_MT2_0_0_MT2_MT2_MU2_1   = qli4_(&zero,&mt2,&zero,&mt2,&T11,&T12,&zero,&zero,&mt2,&mt2,&MUR2,&I);
  I4_0_T12_0_T11_MT2_MT2_0_0_MT2_MT2_MU2_1   = qli4_(&zero,&mt2,&zero,&mt2,&T12,&T11,&zero,&zero,&mt2,&mt2,&MUR2,&I);


#ifdef WITH_NON_FACT_DIAGRAMS
  setlambda(0);

  
  setlambda(-1);

#endif
 }




void Print_SI(int eps)
{
  std::cout << std::endl << " Current values of scalar integrals ";
  switch(eps)
    {
    case 0:
      std::cout << "[fin.]" << std::endl;
      // 1-point functions
      PRINT(I1_MT2_MU2_0);
      // 2-point functions  
      PRINT(I2_MT2_0_MT2_MU2_0);  
      PRINT(I2_S12_0_0_MU2_0);  
      PRINT(I2_S12_MT2_MT2_MU2_0);
      PRINT(I2_T11_0_MT2_MU2_0);
      PRINT(I2_T12_0_MT2_MU2_0);
      // 3-point functions  
      PRINT(I3_0_0_S12_0_0_0_MU2_0);  
      PRINT(I3_0_T11_MT2_0_0_MT2_MU2_0);
      PRINT(I3_0_T12_MT2_0_0_MT2_MU2_0);  
      PRINT(I3_MT2_0_T11_0_MT2_MT2_MU2_0);  
      PRINT(I3_MT2_0_T12_0_MT2_MT2_MU2_0);
      PRINT(I3_MT2_S12_MT2_0_MT2_MT2_MU2_0); 
      PRINT(I3_S12_0_0_MT2_MT2_MT2_MU2_0);  
      PRINT(I3_S12_MT2_MT2_0_0_MT2_MU2_0);  
      // 4-point functions
      PRINT(I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_0);
      PRINT(I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_0);
      break;
    case 1:
      std::cout << "[eps^1 poles]" << std::endl;
      // 1-point functions
      PRINT(I1_MT2_MU2_1);
      // 2-point functions  
      PRINT(I2_MT2_0_MT2_MU2_1);  
      PRINT(I2_S12_0_0_MU2_1);  
      //PRINT(I2_S12_MT2_MT2_MU2_1);
      PRINT(I2_T11_0_MT2_MU2_1);
      PRINT(I2_T12_0_MT2_MU2_1);
      // 3-point functions  
      PRINT(I3_0_0_S12_0_0_0_MU2_1);  
      PRINT(I3_0_T11_MT2_0_0_MT2_MU2_1);
      PRINT(I3_0_T12_MT2_0_0_MT2_MU2_1);  
      //PRINT(I3_MT2_0_T11_0_MT2_MT2_MU2_1);  
      //PRINT(I3_MT2_0_T12_0_MT2_MT2_MU2_1);
      PRINT(I3_MT2_S12_MT2_0_MT2_MT2_MU2_1); 
      //PRINT(I3_S12_0_0_MT2_MT2_MT2_MU2_1);  
      //PRINT(I3_S12_MT2_MT2_0_0_MT2_MU2_1);  
      // 4-point functions
      PRINT(I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_1);
      PRINT(I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_1);      
      break;
    case 2:
      std::cout << "[eps^2 poles]" << std::endl;
      // 1-point functions
      //PRINT(I1_MT2_MU2_2);
      // 2-point functions  
      //PRINT(I2_MT2_0_MT2_MU2_2);  
      //PRINT(I2_S12_0_0_MU2_2);  
      //PRINT(I2_S12_MT2_MT2_MU2_2);
      //PRINT(I2_T11_0_MT2_MU2_2);
      //PRINT(I2_T12_0_MT2_MU2_2);
      // 3-point functions  
      //PRINT(I3_0_0_S12_0_0_0_MU2_2);  
      //PRINT(I3_0_T11_MT2_0_0_MT2_MU2_2);
      //PRINT(I3_0_T12_MT2_0_0_MT2_MU2_2);  
      //PRINT(I3_MT2_0_T11_0_MT2_MT2_MU2_2);  
      //PRINT(I3_MT2_0_T12_0_MT2_MT2_MU2_2);
      //PRINT(I3_MT2_S12_MT2_0_MT2_MT2_MU2_2); 
      //PRINT(I3_S12_0_0_MT2_MT2_MT2_MU2_2);  
      //PRINT(I3_S12_MT2_MT2_0_0_MT2_MU2_2);  
      // 4-point functions
      //PRINT(I4_MT2_0_0_MT2_T11_S12_0_MT2_MT2_MT2_MU2_2);
      //PRINT(I4_MT2_0_0_MT2_T12_S12_0_MT2_MT2_MT2_MU2_2);      
      break;
    default:
      std::cout << std:: endl << "the value eps = " << eps << " is not suppoorted by this function.";
    }
  std::cout << std::endl;
}



