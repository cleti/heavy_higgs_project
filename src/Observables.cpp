

#include "../inc/Observables.h"




DD_spec DDs_unpolarized[10] = 
  {
    {1,&EvalMW,
     "M_(W)",
     "M_{W}",
     " #Gamma^{-1} #frac{d #Gamma}{d M_{W}}",
     0.01, 100,
     0.9 , 1.1, 
     0.42, 0.5},
    {3,&EvalEW,
     "E_(W)",
     "E_{W}",
     " #Gamma^{-1} #frac{d #Gamma}{d E_{W}}",
     0.01, 100,
     0.8 , 1.1, 
     0.5 , 0.62},
    {0,&EvalEl,
     "E_(l)",
     "E_{l}",
     " #Gamma^{-1} #frac{d #Gamma}{d E_{l}}",
     0.01, 20,
     0.8 , 1.1, 
     0.1 , 0.45},
    {4,&EvalMbl,
     "M_(l,b)",
     "M_{l,b}",
     " #Gamma^{-1} #frac{d #Gamma}{d M_{l,b}}",
     0.01, 20,
     0.8 , 1.1, 
     0.0 , 1.0 },
    {0,&EvalEJb,
     "E_(b-jet)",
     "E_{b}",
     "#Gamma)^{-1} #frac{d #Gamma}{d E_{b}}",
     0.01, 20,
     0.8 , 1.2, 
     0.0 , 0.4 },
    {0,&EvalEJ2,
     "E_(2nd-jet)",
     "E_{2}",
     "#Gamma^{-1} #frac{d #Gamma}{d E_{2}}",
     0.01, 20,
     0.8 , 1.2, 
     0.0 , 0.4 },
    {2,&EvalCJbl,
     "CosT_(b-jet,l)",
     "cos(#theta_{bl})",
     "(#Gamma)^{-1} #frac{d #Gamma}{d cos(#theta_{bl})}",
     0.15, 1.05,
     0.97, 1.03, 
     -1.0, 1.0 },
    {2,&EvalCJ2l,
     "CosT_(2nd-jet,l)",
     "cos(#theta_{2l})",
     "#Gamma^{-1} #frac{d #Gamma}{d cos(#theta_{2l})}",
     0.15, 1.05,
     0.97, 1.03, 
     -1.0, 1.0 },
    {2,&EvalHelAngle,
     "CosT_(W,l*)",
     "cos(#theta_{Wl*})",
     "#Gamma^{-1} #frac{d #Gamma}{d cos(#theta_{Wl*})}",
     0.05, 1,
     0.97, 1.03, 
     -1.0, 1.0 },
    {2,&EvalCWb,
    // {1,&EvalMlj,
     "CosT_(W,b)",
     "cos(#theta_{Wb})",
     "(m_{t} #Gamma)^{-1} #frac{d #Gamma}{d cos(#theta_{Wb})}",
     0.02, 8,
     0.85, 1.1, 
     -1.0, 1.0 }
  };

DD_spec DDs_polarized[3] = 
  {
    {2,&EvalCSl,
     "CosT_(S,l)",
     "cos{#theta_{Sl}}",
     "(m_{t} #Gamma)^{-1} #frac{d #Gamma}{d cos{#theta_{Sl}}",
     0.0, 0.0,
     0.0, 0.0, 
    -1.0, 1.0 },
    {2,&EvalCSW,
     "CosT_(S,W)",
     "cos{#theta_{SW}}",
     "(m_{t} #Gamma)^{-1} #frac{d #Gamma}{d cos{#theta_{SW}}",
     0.0, 0.0,
     0.0, 0.0,
    -1.0, 1.0 },
    {2,&EvalCSb,
     "CosT_(S,b)",
     "cos{#theta_{Sb}}",
     "(m_{t} #Gamma)^{-1} #frac{d #Gamma}{d cos{#theta_{Sb}}",
     0.0, 0.0,
     0.0, 0.0, 
    -1.0, 1.0 }
  };   






double EvalMW(C_type C) {
  return sqrt(qW_qW);
}



double EvalEW(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
    case PS4P_B_BB:
      {
	using namespace NSP_PS4P;
	ret = 1.0-kt_k1-kt_k2-kt_k3;
      }
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = 1.0-kt_kE-kt_kS;
      }
      break;
    case PS3P_R_G_13:
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = 1.0-kt_kE-kt_kS;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = 1.0-kt_kE-kt_kS;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = 1.0-kt_kb-kt_kg;
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	ret = 1.0-kt_kb;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalEW,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}


double EvalHelAngle(C_type C) {
  return CosTWl;
}
double EvalCSl(C_type C) {
  return CosTSl;
}
double EvalCSW(C_type C) {
  return CosTSW;
}


// lepton energy
double EvalEl(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
    case PS4P_B_BB:
      {
	using namespace NSP_PS4P; 
	ret = kt_kl;
      }
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = kt_kL;
      }
      break;
    case PS3P_R_G_13:
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13; 
	ret = kt_kL;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = kt_kL;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = kt_kl;
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	ret = kt_kl;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalEl,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}





////////////////////////////////////////////////////////////////////////////////
////////////////////// PARTON JET OBSERVABLES //////////////////////////////////
//////////////////// (reco. scheme dependent!!!) ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// (t-restframe) energy of b flavoured parton jet
double EvalEJb(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:
	  ret = kt_k2+kt_k3;
	  break;
	case 2:
	  ret = kt_k1+kt_k3;
	  break;
	case 3:
	  ret = kt_k3;
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1: // parton 1: b-bar -> (23)= bb-jet => that would be 2 b-flavoured jets
	  ret = 0.0;
	  WARN(EvalEJb,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:
	  ret = kt_k2;
	  break;
	case 3:
	  ret = kt_k3;
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = kt_kS;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = kt_kS;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = kt_kE;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = kt_kE;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = kt_kb;
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	ret = kt_kb;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalEJb,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}


// (t-restframe) energy of unflavoured parton jet
double EvalEJ2(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:
	  ret = kt_k1;
	  break;
	case 2:
	  ret = kt_k2;
	  break;
	case 3:
	  ret = kt_k1+kt_k2;
	  break;
	}
      break;
    case PS4P_B_BB:
      using namespace NSP_PS4P;
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:
	  ret = 0;
	  WARN(EvalEJ2,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:
	  ret = kt_k1+kt_k3;
	  break;
	case 3:
	  ret = kt_k1+kt_k2;
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = kt_kE;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = kt_kE;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = kt_kS;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = kt_kS;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = kt_kg;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalEJ2,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}






// invariant mass lepton + b-jet
double EvalMbl(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
      using namespace NSP_PS4P;
	case 1: // k1 is resolved => Kb = k2+k3
	  ret = sqrt(2.0*(k2_k3+k2_kl+k3_kl));
	  break;
	case 2: // k2 is resolved => Kb = k1+k3
	  ret = sqrt(2.0*(k1_k3+k1_kl+k3_kl));
	  break;
	case 3: // k3 is resolved => Kb = k3
	  ret = sqrt(2.0*(k3_kl));
	  break;
	}
      break;
    case PS4P_B_BB:
      using namespace NSP_PS4P;
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1: // k1 resolved => 2 b-jets!!!
	  ret = 0;
	  WARN(EvalMbl,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2: // k2 resolved => Kb = k2
	  ret = sqrt(2.0*(k2_kl));
	  break;
	case 3: // k3 resolved => Kb = k3
	  ret = sqrt(2.0*(k3_kl));
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = sqrt(2.0*(kS_kL));
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = sqrt(2.0*(kS_kL));
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = sqrt(2.0*(kE_kL));
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = sqrt(2.0*(kE_kL));
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = sqrt(2.0*(kb_kl));
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	ret = sqrt(2.0*(kb_kl));
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalMbl,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}



// invariant mass lepton + 2nd-jet
double EvalM2l(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
      using namespace NSP_PS4P;
	case 1: // k1 is resolved => K2 = k1
	  ret = sqrt(2.0*(k1_kl));
	  break;
	case 2: // k2 is resolved => K2 = k2
	  ret = sqrt(2.0*(k2_kl));
	  break;
	case 3: // k3 is resolved => K2 = k1+k2
	  ret = sqrt(2.0*(k1_kl+k2_kl+k1_k2));
	  break;
	}
      break;
    case PS4P_B_BB:
      using namespace NSP_PS4P;
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1: // k1 resolved => 2 b-jets!!!
	  ret = 0;
	  WARN(EvalMbl,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2: // k2 resolved => K2 = k1+k3
	  ret = sqrt(2.0*(k1_kl+k3_kl+k1_k3));
	  break;
	case 3: // k3 resolved => K2 = k1+k2
	  ret = sqrt(2.0*(k1_kl+k2_kl+k1_k2));
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = sqrt(2.0*(kE_kL));
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = sqrt(2.0*(kE_kL));
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = sqrt(2.0*(kS_kL));
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = sqrt(2.0*(kS_kL));
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = sqrt(2.0*(kg_kl));
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalM2l,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}








// angle of parton b-jet and lepton
double EvalCJbl(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:// b-jet = partons 2 and 3
	case 2:// b-jet = partons 1 and 3
	  ret = ( (kt_kl*kt_Kij - Kij_kl) / (kt_kl*sqrt(pow(kt_Kij,2.0)-Kij_Kij)) );
	  break;
	case 3:// b-jet = parton 3
	  ret = ( (kt_kl*kt_Kk - Kk_kl) / (kt_kl*kt_Kk) );
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:
	  ret = 0;
	  WARN(EvalCJbl,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// b-jet = parton 2
	case 3:// b-jet = parton 3
	  ret = ( (kt_kl*kt_Kk - Kk_kl) / (kt_kl*kt_Kk) );
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = 1.0-kS_kL/kt_kL/kt_kS;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = 1.0-kS_kL/kt_kL/kt_kS;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = 1.0-kE_kL/kt_kL/kt_kE;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = 1.0-kE_kL/kt_kL/kt_kE;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = 1.0-CosTbl;
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	ret = 1.0-CosTbl;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalCJbl,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}

// angle of parton non-b-jet and lepton
double EvalCJ2l(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:// g-jet = partons 1
	case 2:// g-jet = partons 2
	  ret = ( (kt_kl*kt_Kk - Kk_kl) / (kt_kl*kt_Kk) );
	  break;
	case 3:// g-jet = parton 1+2
	  ret = ( (kt_kl*kt_Kij - Kij_kl) / (kt_kl*sqrt(pow(kt_Kij,2.0)-Kij_Kij)) );
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	  using namespace NSP_PS4P;
	case 1:
	  ret = 0;
	  WARN(EvalCJ2l,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// g-jet = parton 1+3
	case 3:// g-jet = parton 1+2
	  ret = ( (kt_kl*kt_Kij - Kij_kl) / (kt_kl*sqrt(pow(kt_Kij,2.0)-Kij_Kij)) );
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = 1.0-kE_kL/kt_kL/kt_kE;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = 1.0-kE_kL/kt_kL/kt_kE;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = 1.0-kS_kL/kt_kL/kt_kS;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = 1.0-kS_kL/kt_kL/kt_kS;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = 1.0-CosTgl;}
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalCJ2l,"C_type not specified");
	  MSG = true;
	}
    }
  return ret;
}



// W+, b-jet angle
double EvalCWb(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	case 1:// b-jet = partons 2 and 3
	case 2:// b-jet = partons 1 and 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kij_qW = kt_Kij - Kij_Kk - Kij_Kij;
	    ret = (kt_qW*kt_Kij - Kij_qW) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij)*(pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kk_qW = kt_Kk - Kij_Kk;
	    ret = (kt_qW*kt_Kk  - Kk_qW ) / kt_Kk / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	case 1:
	  ret = 0;
	  WARN(EvalCWb,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// b-jet = parton 2
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kk_qW = kt_Kk - Kij_Kk;
	    ret = (kt_qW*kt_Kk  - Kk_qW ) / kt_Kk / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kS_qW = kt_kS - kE_kS;
	ret = (kt_qW*kt_kS  - kS_qW ) / kt_kS / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kS_qW = kt_kS - kE_kS;
	ret = (kt_qW*kt_kS  - kS_qW ) / kt_kS / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kE_qW = kt_kE - kE_kS;
	ret = (kt_qW*kt_kE  - kE_qW ) / kt_kE / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kE_qW = kt_kE - kE_kS;
	ret = (kt_qW*kt_kE  - kE_qW ) / kt_kE / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	double kt_qW = 1.0-kt_kb-kt_kg;
	double kb_qW = kt_kb - kg_kb;
	ret = (kt_qW*kt_kb  - kb_qW ) / kt_kb / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS2P:
      {
	using namespace NSP_PS2P;
	double kt_qW = 1.0-kt_kb;
	ret = (kt_qW*kt_kb  - kt_kb ) / kt_kb / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    default:
      if (!MSG) {
	WARN(EvalCWb,"C_type not specified");
	MSG = true;
      }
    }
  return ret;
}


// W+, 2nd-jet angle
double EvalCW2(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	case 1:// b-jet = partons 2 and 3
	case 2:// b-jet = partons 1 and 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kk_qW = kt_Kk - Kij_Kk;
	    ret = (kt_qW*kt_Kk  - Kk_qW ) / kt_Kk / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kij_qW = kt_Kij - Kij_Kk - Kij_Kij;
	    ret = (kt_qW*kt_Kij - Kij_qW) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij)*(pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	case 1:
	  ret = 0;
	  WARN(EvalCW2,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// b-jet = parton 2
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    double kt_qW = 1.0-kt_k1-kt_k2-kt_k3;
	    double Kij_qW = kt_Kij - Kij_Kk - Kij_Kij;
	    ret = (kt_qW*kt_Kij - Kij_qW) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij)*(pow(kt_qW,2.0)-qW_qW) ) ;
	  }
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kE_qW = kt_kE - kE_kS;
	ret = (kt_qW*kt_kE  - kE_qW ) / kt_kE / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kE_qW = kt_kE - kE_kS;
	ret = (kt_qW*kt_kE  - kE_qW ) / kt_kE / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kS_qW = kt_kS - kE_kS;
	ret = (kt_qW*kt_kS  - kS_qW ) / kt_kS / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	double kt_qW = 1.0-kt_kS-kt_kE;
	double kS_qW = kt_kS - kE_kS;
	ret = (kt_qW*kt_kS  - kS_qW ) / kt_kS / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	double kt_qW = 1.0-kt_kb-kt_kg;
	double kg_qW = kt_kg - kg_kb;
	ret = (kt_qW*kt_kg  - kg_qW ) / kt_kg / sqrt( (pow(kt_qW,2.0)-qW_qW) ) ;
      }
      break;
    default:
      if (!MSG) {
	WARN(EvalCW2,"C_type not specified");
	MSG = true;
      }
    }
  return ret;
}



// PROBLEM>????
// angle between top spin and b-jet direction
double EvalCSb(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	case 1:// b-jet = partons 2 and 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k2-S_k3 ) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij) ) ;
	  }
	  break;
	case 2:// b-jet = partons 1 and 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k1-S_k3 ) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij) ) ;
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k3 ) / kt_k3;
	  }
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	case 1:
	  ret = 0;
	  WARN(EvalCSb,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// b-jet = parton 2
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k2 ) / kt_k2;
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k3 ) / kt_k3;
	  }
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = ( - S_kS ) / kt_kS;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = ( - S_kS ) / kt_kS;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = ( - S_kE ) / kt_kE;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = ( - S_kE ) / kt_kE;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = ( - S_kb ) / kt_kb;
      }
      break;
    default:
      if (!MSG) 
	{
	  WARN(EvalCSb,"C_type not specified");
	  MSG = true;
	}
    }

  CHECKANGLE(1.0-ret);
  return ret;
}


// angle between top spin and 2nd-jet direction
double EvalCS2(C_type C) {
  static bool MSG = false;
  double ret;
  switch (C)
    {
    case PS4P_B_GG:
    case PS4P_B_QQ:
      switch (RESOLVED_PARTON)
	{
	case 1:// b-jet = partons 2 and 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k1 ) / kt_k1;
	    if (ret<=-1. || ret>= 1.) WARN(EvalCS2,"angle out of range");
	  }
	  break;
	case 2:// b-jet = partons 1 and 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k2 ) / kt_k2;
	    if (ret<=-1. || ret>= 1.) WARN(EvalCS2,"angle out of range");
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k1-S_k2 ) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij) ) ;
	    if (ret<=-1. || ret>= 1.) WARN(EvalCS2,"angle out of range");
	  }
	  break;
	}
      break;
    case PS4P_B_BB:
      switch (RESOLVED_PARTON)
	{
	case 1:
	  ret = 0;
	  WARN(EvalCS2,"This configuration corresponds to 2 b-flavoured jets !!!");
	  break;
	case 2:// b-jet = parton 2
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k1-S_k3 ) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij) ) ;
	    if (ret<=-1. || ret>= 1.) WARN(EvalCS2,"angle out of range");
	  }
	  break;
	case 3:// b-jet = parton 3
	  {
	    using namespace NSP_PS4P;
	    ret = (-S_k1-S_k2 ) / sqrt( (pow(kt_Kij,2.0)-Kij_Kij) ) ;
	    if (ret<=-1. || ret>= 1.) WARN(EvalCS2,"angle out of range");
	  }
	  break;
	}
      break;
    case PS3P_R_G_12:
      {
	using namespace NSP_PS3P_R_12;
	ret = ( - S_kE ) / kt_kE;
      }
      break;
    case PS3P_R_G_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = ( - S_kE ) / kt_kE;
      }
      break;
    case PS3P_R_B_13:
      {
	using namespace NSP_PS3P_R_13;
	ret = ( - S_kS ) / kt_kS;
      }
      break;
    case PS3P_R_B_23:
      {
	using namespace NSP_PS3P_R_23;
	ret = ( - S_kS ) / kt_kS;
      }
      break;
    case PS3P_B_G:
      {
	using namespace NSP_PS3P;
	ret = ( - S_kg ) / kt_kg;
      }
      break;
    default:
      if (!MSG) {
	WARN(EvalCS2,"C_type not specified");
	MSG = true;
      }
    }
  return ret;
}





////////////////////////////////////////////////////////////////////////////////
////////////////////// PARTON JET OBSERVABLES //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
