


#include <iostream>
#include <cmath>
#include <bitset>
#include <random>
#include <memory>
using namespace std;

#include "../inc/Functions_Shared.h"
#include "../inc/Makros.h"

#include "../inc/Functions_pp_ttX_V.h"
#include "../inc/Functions_pp_ttX_ID.h"
#include "../inc/Functions_pp_ttX_R.h"
#include "../inc/Functions_pp_ttX_UID.h"
#include "../inc/Integrands_pp_ttX.h"
#include "../inc/Integrator.h"
#include "../inc/HiggsModel.h"
#include "../inc/ScalarIntegrals.h"

#include "../ext/LoopTools-2.12/build/clooptools.h"

#ifdef __cplusplus
extern "C" {
#endif
  extern void qlinit_();
  extern c_double qli1c_(c_double*,double*,int*);
#ifdef __cplusplus
}
#endif



double fnc1(double* x, size_t dim, void* par)
{
  return sqrt(x[0]);
}

double fnc2(double* x, size_t dim, void* par)
{
  return x[1]+0.5*sqrt(x[0]);
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  using namespace Constants;
  using RunParameters::g_flags_eval_v;
  using RunParameters::g_flags_eval_id;
  using RunParameters::g_flags_eval_r; 

  qlinit_();
  ltini();

  int PREC = (argc>1) ? atoi(argv[1]) : 5;
  cout << setprecision(PREC);
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // input parameters ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  double MT = 173.5;
  const double S = 5.0;
  int GGH_EFF = 0;
  
  HiggsModel THDM;
  // from now on all dimensionfull quantities are normalized to the top mass
  THDM.SetScale(MT);
  THDM.SetAlphaS(0.105117798954049);
  THDM.SetMUR(2.0*MT);
  THDM.SetMUF(2.0*MT);
  THDM.SetVH(246.0); 
  THDM.SetMt(MT);
  THDM.SetMb(0.0);

  
  THDM.AddBoson(
  		500.0,
  		3.245439053359254e+01,
  		1,1,
  		0,0);
  // THDM.AddBoson(
  // 		600.0,
  // 		67.71,
  // 		-1.19,1.18,
  // 		0,0);

  THDM.SetHiggsPrefactors(S,GGH_EFF);
  THDM.Print(cout,MT);
  double const& mt2 = THDM.mt2();
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // {
  //   using namespace HiggsBosons;
  //   SetPhi1(500.0,
  //   	    3.245439053359254e+01,
  //   	    1.0,
  //   	    1.0,
  //   	    0.0,
  //   	    0.0
  //   	    );
  //   SetPhi2(600.0,
  //   	    67.71,
  //   	    -1.19,
  //   	    1.18,
  //   	    0.0,
  //   	    0.0);
  //   RunParameters::TwoHDM = 1;
  //   ResetHiggsPrefactors(S,GGH_EFF);
  //   PrintHiggsPrefactors();
  // }
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // born + virtual + dipoles  integration
  {
    PS_2_2  ps_gg_tt(mt2,mt2,"gg->tt-bar");
    PS_2_2  ps_gg_tt_x(mt2,mt2,"gg->tt-bar (x)");
    // set flags to specify which matrix elements will be evaluated
    USET_EVAL_ALL(g_flags_eval_v);
    // Born
    SET_EVAL_B_PHIxPHI(g_flags_eval_v);
    USET_EVAL_B_PHIxQCD(g_flags_eval_v);
    USET_EVAL_B_QCDxQCD(g_flags_eval_v);
    // virtual corrections
    USET_EVAL_V_PHIxPHI(g_flags_eval_v);
    SET_EVAL_V_PHI1xQCD0(g_flags_eval_v);
    USET_EVAL_V_PHI0xQCD1(g_flags_eval_v);
    //SET_EVAL_V_NF(g_flags_eval_v);

    
    // integrated dipoles
    SET_EVAL_B_PHIxQCD(g_flags_eval_id);
    SET_EVAL_B_PHIxPHI(g_flags_eval_id);
    // integrated dipoles 
    SET_EVAL_D_DELTA(g_flags_eval_id);
    SET_EVAL_D_CONT(g_flags_eval_id);
    SET_EVAL_D_END(g_flags_eval_id);

    cout << endl << " EVAL_V_FLAGS = " << bitset<16>(g_flags_eval_v ).to_string() << endl;
    cout << endl << " EVAL_D_FLAGS = " << bitset<16>(g_flags_eval_id).to_string() << endl;

    
    double res_b=0,res_v=0,res_d=0,res_dx=0;
    double y_values[] = {-0.9};//-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9};
    double x_values[] = {0.85};//,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};

    for (auto x : x_values)
      {
	for (auto y : y_values)
	  {
	    if (ps_gg_tt.set(sqrt(S),y) )
	      {

#ifdef WITH_T_SPIN
		cout << endl << "Using spin dependent matrix elements. Spin 4-vectors are:" << endl;
		// ALWAYS use spin vectors with the properties
		// s_i.s_i = -1 and s_i.k_i = 0 !!!
		// otherwise tests will fail, since these eqs. are used in the MEs
		LT boost;
		FV& s1 = ps_gg_tt.s1();
		FV& s2 = ps_gg_tt.s2();
		//////////////////////////////////////
		s1 = {0.0,0.0,sqrt(0.5),+sqrt(0.5)};
		s2 = {0.0,sqrt(0.5),0.0,-sqrt(0.5)};
		//////////////////////////////////////
		FV const& k1 = ps_gg_tt.k1();
		FV const& k2 = ps_gg_tt.k2();
		boost.set_boost(k1,1);
		boost.apply(s1);
		boost.set_boost(k2,1);
		boost.apply(s2);
		PRINT_4VEC(s1);
		PRINT_4VEC(s2);
		PRINT(sp(s1,k1));
		PRINT(sp(s2,k2));
#endif

		
		ps_gg_tt.set_x(x);
		ps_gg_tt_x = ps_gg_tt;
		ps_gg_tt_x.set(sqrt(x*S),y);
		cout << endl;
		// cout << endl << " Setting '" << ps_gg_tt.d_name << "' to ";
		cout << endl << "sqrt(s) = " << ps_gg_tt.get_rs() << "*mScale , y = " << y << ", x = " << x;
		cout << endl << "---------------------------------------------------";
		cout << endl << " Born                   = " << (res_b = PREF_GG*Eval_B   (ps_gg_tt,THDM,g_flags_eval_v,GGH_EFF)) << "  " << endl;
		cout << endl << " Virt.                  = " << (res_v = PREF_GG*Eval_V   (ps_gg_tt,THDM,g_flags_eval_v)) << "  " << endl;
		cout << endl << " DIP                    = " << (res_d = PREF_GG*Eval_ID  (ps_gg_tt,THDM,g_flags_eval_id)) << "  " << endl;		
		cout << endl << " DIP_X                  = " << (res_dx= PREF_GG*Eval_ID_X(ps_gg_tt_x,THDM,g_flags_eval_id)) << "  " << endl;
		cout << endl << " DIP Sum                = " << res_d+res_dx;
		cout << endl << "---------------------------------------------------";
		cout << endl << " NLO                    = " << (res_d+res_dx+res_v+res_b) << "  " << endl;
	      }
	  }
	cout << endl << " =========================================================================================" << endl;
      }
  }

  
//   // exit(1);
  
  // test real corrections
  {
    PS_2_3 ps_gg_ttg(mt2,mt2,0.0,"gg->tt-bar+g");

    double res_r_phi  = 0.0;
    double res_r_int  = 0.0;
    double res_d_phi  = 0.0;
    double res_d_int  = 0.0;
    double sum_r = 0.0;
    double sum_d = 0.0;
    double rs_part= 5.0;
    double y_cmf  = 0.5;// 0.5!!!
    double phi_cmf= 0.0;
    double M_tt   = 3.0;
    double y_tt   = 0.0;
    double phi_tt = 0.0;
    FV n;
  
    if ( ps_gg_ttg.set(rs_part,y_cmf,phi_cmf,M_tt,y_tt,phi_tt) )
      {
	double const& mScale2 = THDM.Scale2();
	
	// FV& p1 = ps_gg_ttg.p1();
	// FV& p2 = ps_gg_ttg.p2();
	// FV& p3 = ps_gg_ttg.p3();
	// n = {p3[0],-p3[1],-p3[2],-p3[3]};
      
	// double k1k2 = sp(k1,k2);
	// double p1p2 = sp(p1,p2);
	// double k1p1 = sp(k1,p1);
	// double k1p2 = sp(k1,p2);
	// double k2p1 = sp(k2,p1);
	// double k2p2 = sp(k2,p2);
	// double p1q  = sp(p1,p3);
	// double p2q  = sp(p2,p3);
	// double k1q  = sp(k1,p3);
	// double k2q  = sp(k2,p3);
	// double k1n  = sp(k1,n);
	// double k2n  = sp(k2,n);
	// double p1n  = sp(p1,n);
	// double p2n  = sp(p2,n);
	// double qn   = sp(p3,n);
			
	// cout << setprecision(15);
      
	// PRINT(k1k2);
	// PRINT(p1p2);
	// PRINT(k1p1);
	// PRINT(k1p2);
	// PRINT(k2p1);
	// PRINT(k2p2);
	// PRINT(p1q);
	// PRINT(p2q);
	// PRINT(k1q);
	// PRINT(k2q);
	// PRINT(k1n);
	// PRINT(k2n);
	// PRINT(p1n);
	// PRINT(p2n);
	// PRINT(qn);
		
	// ps_gg_ttg.print();
	// interference terms
	USET_EVAL_R_ALL(g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_ISR(g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_FSR(g_flags_eval_r);
	SET_EVAL_R_ISR_ISR(g_flags_eval_r);
	// USET_EVAL_R_ISR_ISR(g_flags_eval_r);
	// USET_EVAL_R_FSR_ISR(g_flags_eval_r);
	// USET_EVAL_R_FSR_INT(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r ).to_string() << endl;
	cout << endl << " 2*RE[M_QCD * M_phi ]     = " << (res_r_int = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip 2*RE[M_QCD * M_phi ] = " << (res_d_int = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	// phi^2 terms [FSR]
	USET_EVAL_R_ALL(g_flags_eval_r);
	SET_EVAL_R_PHIxPHI_FSR(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r ).to_string() << " [FSR] " << endl;
	cout << endl << " [M_phi]^2            = " << (res_r_phi = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << (res_d_phi = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	sum_r += res_r_phi;
	sum_d += res_d_phi;
	// phi^2 terms [ISR]
	USET_EVAL_R_ALL(g_flags_eval_r);
	SET_EVAL_R_PHIxPHI_ISR(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r ).to_string() << " [ISR] " << endl;
	cout << endl << " [M_phi]^2            = " << (res_r_phi = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << (res_d_phi = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	sum_r += res_r_phi;
	sum_d += res_d_phi;

	// qq
	USET_EVAL_R_ALL(g_flags_eval_r);
	SET_EVAL_R_PHIxPHI_QQ(g_flags_eval_r);
	USET_EVAL_R_PHIxQCD_QQ(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r ).to_string() << " [QQ] " << endl;
	cout << endl << " [M_phi]^2            = " << (res_r_phi = PREF_QQ*Eval_R_QQ(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;

	// qg
	USET_EVAL_R_ALL(g_flags_eval_r);
	SET_EVAL_R_PHIxPHI_QG(g_flags_eval_r);
	USET_EVAL_R_PHIxQCD_QG(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<16>(g_flags_eval_r ).to_string() << " [QG] " << endl;
	cout << endl << " [M_phi]^2            = " << (res_r_phi = PREF_QG*Eval_R_QG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << (res_d_phi = PREF_QG*Eval_UID_QG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;

	cout << endl << " === SWAP !!! === " << endl;
        ps_gg_ttg.swap();
	cout << endl << " [M_phi]^2            = " << (res_r_phi = PREF_QG*Eval_R_QG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << (res_d_phi = PREF_QG*Eval_UID_QG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	ps_gg_ttg.swap();

	cout << endl << " sum_r = " << sum_r/mScale2;
	cout << endl << " sum_d = " << sum_d/mScale2 << endl;
      }
  }
  return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
