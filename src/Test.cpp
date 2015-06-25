


#include <iostream>
#include <cmath>
#include <bitset>
#include <random>
#include <memory>

#include <sys/ioctl.h>
#include <cstdio>
//#include <unistd.h>

using namespace std;

#include "../inc/Functions_Shared.h"
#include "../inc/Makros.h"

#include "../inc/Functions_pp_ttX_V.h"
#include "../inc/Functions_pp_ttX_ID.h"
#include "../inc/Functions_pp_ttX_R.h"
#include "../inc/Functions_pp_ttX_UID.h"

#ifdef WITH_T_SPIN
#include "../inc/Functions_tDecay.h"
#endif

#include "../inc/Integrands_pp_ttX.h"
#include "../inc/Integrator.h"
#include "../inc/HiggsModel.h"
#include "../inc/ScalarIntegrals.h"

#include "../ext/LoopTools-2.12/build/clooptools.h"

// LHAPDF header files
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/Info.h"
#include "LHAPDF/Config.h"
#include "LHAPDF/PDFInfo.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Factories.h"

#include <boost/type_traits.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  extern void qlinit_();
  extern c_double qli1c_(c_double*,double*,int*);
#ifdef __cplusplus
}
#endif


ulong g_flags_eval_v;
ulong g_flags_eval_id;
ulong g_flags_eval_r;


double fnc1(double* x, size_t dim, void* par)
{
  return sqrt(x[0]);
}

double fnc2(double* x, size_t dim, void* par)
{
  return x[1]+0.5*sqrt(x[0]);
}


static double Ks(double const& as,
		 double const& mu,
		 double const& mh,
		 double const& mt,
		 double const& nf = 6)
{
  using namespace Constants;
  return 1.0+as/Pi*(95.0/4.0-7.0/6.0*nf+(33.0-2.0*nf)/3.0*log(mu/mh))+pow(as/Pi,2)*(156.808-2.0*5.708*log(mt/mh));
}
static double Kp(double const& as,
		 double const& mu,
		 double const& mh,
		 double const& mt,
		 double const& nf = 6)
{
  using namespace Constants;
  return 1.0+as/Pi*(97.0/4.0-7.0/6.0*nf+(33.0-2.0*nf)/3.0*log(mu/mh))+pow(as/Pi,2)*(171.544-2.0*5.*log(mt/mh));
}

// computes the phi -> gg decay rate inlcuding NNLO QCD corrections in the large mt limit if NNLO = true
static double Gamma_hgg(
			int h,
			HiggsModel& model,
			double const& mt_MS,			
			bool NNLO
			)
{
  using namespace Constants;

  double const& as     = model.AlphaS();
  double const& mt_pole= model.mt();
  double const& mb_MS  = model.mb();
  HPtr phi = model.GetBoson(h);
  double const& mu     = model.MUR();//phi->M();

    PRINT(mu);
    PRINT(as);
    PRINT(mt_pole);
    PRINT(mt_MS);
    PRINT(mb_MS);
  
  double ks = 1.0;
  double kp = 1.0;
  if (NNLO)
    {
      int nf = 5;
      // running top mass for the radiative corrections
      ks = Ks(as,mu,phi->M(),mt_MS,nf);
      kp = Kp(as,mu,phi->M(),mt_MS,nf);
    }

  // top pole mass and bottom MS-bar mass  for evaluation of the 1-loop form factors
  phi->SetFormFactors(phi->M2(),pow(mt_pole,2),pow(mb_MS,2));
  c_double Fs1 = phi->GetFs()*(-4.0)*phi->M2();
  c_double Fp1 = phi->GetFp()*8.0*phi->M2();
  
  return 1.0/(32.0*Pi)*CF*CA/(phi->M())*pow(as/Pi,2)*(ks*norm(Fs1)+kp*norm(Fp1));
}

static double x (const FV& p_i, const FV& p_a, const FV& p_b)
{
  double t1 = sp(p_a, p_b);
  double t2 = sp(p_i, p_a);
  double t3 = sp(p_i, p_b);
  return ((t1 - t2 - t3) / t1);
}





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{  
  using namespace Constants;


  qlinit_();
  //  ltini();

  int PREC = (argc>1) ? atoi(argv[1]) : 5;
  cout << setprecision(PREC);

  typedef std::vector<double> DVec;
  DVec::iterator it;
  boost::is_pointer<DVec::iterator>::type it_type;
  boost::is_pointer<DVec*>::type p_type;
  cout << endl << it_type << endl;
  cout << endl << p_type << endl;

  FV a = {1,2,3,4};
  
  cout << endl << sizeof(a)/sizeof(int) << endl;    
  cout << endl << sizeof(PS_2_1)/sizeof(int) << endl;  
  cout << endl << sizeof(PS_2_2)/sizeof(int) << endl;
  cout << endl << sizeof(PS_2_3)/sizeof(int) << endl;  

  EXIT(1);
  
  {
    int cols = 80;
    int lines = 24;

#ifdef TIOCGSIZE
    struct ttysize ts;
    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
    cols = ts.ts_cols;
    lines = ts.ts_lines;
#elif defined(TIOCGWINSZ)
    struct winsize ts;
    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
    cols = ts.ws_col;
    lines = ts.ws_row;
#endif /* TIOCGSIZE */

    cout << endl << cols << " x " << lines << endl;
    //exit(1);
  }
  

  ////////////////////////////////////////////////////////////////////////////////////////////
  // PDF scan ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  {
    double MUF2 = pow(265.0,2);
    // create the PDFs
    LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10nlo", 0);

    static auto quark_pids     = { 5,  4,  3,  2,  1};
    double x1 = 0.0;
    double x2 = 0.0;
    double start= 0.1;
    double end  = 0.9;
    double step = 0.1;
    int steps = (end-start)/step;

    cout << endl;
    cout << " QQ-bar PDFs Q=[u,d,s,c,b]: ";
    cout << endl;
    cout << setw(10) << "" << "|| ";
    for (x1=start;x1<end;x1+=step)
      {
	cout << setw(10) << setprecision(2) << x1 << " ";
      }
    cout << endl;
    cout << setw(11*(steps+2)+1) << setfill('=') << "";cout << setfill(' ');
    cout << endl;
    for (x1=start;x1<end;x1+=step)
      {
	cout << setw(10) << setprecision(3) << x1 << "|| ";
	for (x2=start;x2<end;x2+=step)
	  {
	    double qq_pdf = 0.0;
	    // get the quark/antiquark pdfs
	    for (auto i: quark_pids)
	      {
		double f1_q  = pdf->xfxQ2(+i, x1, MUF2) / x1;
		double f2_qb = pdf->xfxQ2(-i, x2, MUF2) / x2;
		double f1_qb = pdf->xfxQ2(-i, x1, MUF2) / x1;
		double f2_q  = pdf->xfxQ2(+i, x2, MUF2) / x2;
		qq_pdf += (f1_q*f2_qb+f1_qb*f2_q);
	      }
	    cout << setw(10) << setprecision(3) << qq_pdf << " ";
	  }
	cout << endl;
      }
    cout << endl;

  
    cout << endl;
    cout << " GG PDFs: ";
    cout << endl;
    cout << setw(10) << "" << "|| ";
    for (x1=start;x1<end;x1+=step)
      {
	cout << setw(10) << setprecision(2) << x1 << " ";
      }
    cout << endl;
    cout << setw(11*(steps+2)+1) << setfill('=') << "";cout << setfill(' ');
    cout << endl;
    for (x1=start;x1<end;x1+=step)
      {
	cout << setw(10) << setprecision(3) << x1 << "|| ";
	for (x2=start;x2<end;x2+=step)
	  {
	    double gg_pdf = 1.0;
	    gg_pdf *= pdf->xfxQ2(21, x1, MUF2) / x1;
	    gg_pdf *= pdf->xfxQ2(21, x2, MUF2) / x2;
	    cout << setw(10) << setprecision(3) << gg_pdf << " ";
	  }
	cout << endl;
      }
    cout << endl;

    
    delete pdf;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////


  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // compute H->gg decay widths //////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  {
    cout << setprecision(PREC);

    // create the PDFs
    LHAPDF::PDF* pdf = LHAPDF::mkPDF("CT10nlo", 0);
    
    // running MS-bar parameters at mu= {m2 , m3 , (m2+m3)/2} [polemasses for threshold]
    // double mt[3] = {151.309,152.148,151.719};
    // double mb[3] = {2.62942,2.64399,2.63654};
    // double as[3] = {0.0935236,0.0943685,0.0939362};

    // pole masses
    double mtp = 173.34;
    double mbp = 4.75;

    const int NCOL = 5;
    const int NSCEN= 3;

    // took the values of the MS-bar masses from 2HDMC
    // running MS-bar top mass for the scales relevant in the scenarios 1-3:
    double mt[NCOL*NSCEN] = {
      //mu=m_2,  m3,        (m2+m3)/4,   (m2+m3)/2,   (m2+m3)/8,   
      151.309,   152.148,   160.010,     151.719,     169.347,  // scenario 1
      151.309,   148.715,   157.923,     149.918,     167.120,  // scenario 2
      152.37 ,   147.331,   157.438,     149.499,     166.603   // scenario 3 neu     
    };
    // running MS-bar bottom mass for the scales relevant in the scenarios 1-3:
    double mb[NCOL*NSCEN] = {
      //mu=m_2,  m3,        (m2+m3)/4,   (m2+m3)/2,   (m2+m3)/8,   
      2.6294,    2.6440,    2.7806,      2.6365,      2.9429,   // scenario 1
      2.6294,    2.5843,    2.7444,      2.6053,      2.9042,   // scenario 2
      2.6479,    2.5603,    2.7359,      2.5980,      2.8952    // scenario 3 neu     
    };    

  
    HiggsModel THDM("scenario 1");
    THDM.SetVH(246.0);
    // top quark pole mass in all BOrn contributions
    THDM.SetMt(mtp);
    
    int scenario = 2; int &S = scenario;

    PRINT(S);
    switch (S)
      {
      case 0:
	// // scenario 1
	THDM.AddBoson(
		      550.0,
		      0,//not relevant here
		      1.4286,0,
		      -0.7,0);
	THDM.AddBoson(
		      510.0,
		      0,//not relevant here
		      0,1.4286,
		      0,0.7);
	break;
      case 1:
	// scenario 2
	THDM.AddBoson(
		      550.0,
		      0,//not relevant here
		      1.4286,0,
		      -0.7,0);
	THDM.AddBoson(
		      700.0,
		      0,//not relevant here
		      0,1.4286,
		      0,0.7);
	break;
      case 2:
	// scenario 3
	THDM.AddBoson(
		      500.0,
		      0,//not relevant here
		      0.8631,0.9881,
		      -0.6420,0.4842);
	THDM.AddBoson(
		      800.0,
		      0,//not relevant here
		      -1.1572,0.9881,
		      0.3480,0.4842); 
        break;
      }
    
    HPtr phi1 = THDM.GetBoson(0);
    HPtr phi2 = THDM.GetBoson(1);

    double MU = (phi1->M()+phi2->M())/4.0;
    THDM.SetMUR(MU);

    THDM.SetMb(mb[NCOL*S+2]);
    THDM.SetMUR(MU);
    THDM.SetAlphaS(pdf->alphasQ(THDM.MUR()));
    double Gamma_2_gg_NLO = Gamma_hgg(0,THDM,mt[NCOL*S+2],1);
    double Gamma_3_gg_NLO = Gamma_hgg(1,THDM,mt[NCOL*S+2],1);

    
    THDM.SetMb(mb[NCOL*S+3]);
    THDM.SetMUR(MU*2.0);
    THDM.SetAlphaS(pdf->alphasQ(THDM.MUR()));
    double Gamma_2_gg_NLO_h = Gamma_hgg(0,THDM,mt[NCOL*S+3],1);
    double Gamma_3_gg_NLO_h = Gamma_hgg(1,THDM,mt[NCOL*S+3],1);

    
    THDM.SetMb(mb[NCOL*S+4]);
    THDM.SetMUR(MU/2.0);
    THDM.SetAlphaS(pdf->alphasQ(THDM.MUR()));
    double Gamma_2_gg_NLO_l = Gamma_hgg(0,THDM,mt[NCOL*S+4],1);
    double Gamma_3_gg_NLO_l = Gamma_hgg(1,THDM,mt[NCOL*S+4],1);


    PRINT(Gamma_2_gg_NLO);
    PRINT(Gamma_3_gg_NLO);
    PRINT(Gamma_2_gg_NLO_h-Gamma_2_gg_NLO);
    PRINT(Gamma_3_gg_NLO_h-Gamma_3_gg_NLO);   
    PRINT(Gamma_2_gg_NLO_l-Gamma_2_gg_NLO);
    PRINT(Gamma_3_gg_NLO_l-Gamma_3_gg_NLO);


   delete pdf;

   //exit(1);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////


  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // input parameters ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  double MT = 173.5;
  const double S = 25.0;
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

  // 2HDM vergleich mit Peter G.
  THDM.AddBoson(
		500.0,
		32.4543905335925,//not relevant here
		1,1,
		0,0);
  THDM.AddBoson(
  		600.0,
  		67.71,//not relevant here
  		-1.19,1.18,
  		0,0);  
  
  // // scenario 1
  // THDM.AddBoson(
  // 		550.0,
  // 		32.77,//not relevant here
  // 		1.4286,0,
  // 		-0.0,0);
  // THDM.AddBoson(
  // 		510.0,
  // 		48.22,//not relevant here
  // 		0,1.4286,
  // 		0,0.0);
  
  THDM.SetHiggsPrefactors(S,GGH_EFF);
  THDM.Print(cout,MT);
  double const& mt2 = THDM.mt2();
  ////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////
  // born + virtual + dipoles  integration
  {
    PS_2_2  ps_gg_tt(mt2,mt2,"gg->tt-bar");
    PS_2_2  ps_gg_tt_x(mt2,mt2,"gg->tt-bar (x)");
    // set flags to specify which matrix elements will be evaluated
    g_flags_eval_v = 0;
    // Born
    SET_EVAL_B_PHIxPHI(g_flags_eval_v);
    USET_EVAL_B_PHIxQCD(g_flags_eval_v);
    USET_EVAL_B_QCDxQCD(g_flags_eval_v);
    // virtual corrections
    SET_EVAL_V_PHIxPHI(g_flags_eval_v);
    SET_EVAL_V_PHI1xQCD0(g_flags_eval_v);
    SET_EVAL_V_PHI0xQCD1(g_flags_eval_v);
    //SET_EVAL_V_NF(g_flags_eval_v);

    
    // integrated dipoles
    g_flags_eval_id = 0;
    USET_EVAL_B_PHIxQCD(g_flags_eval_id);
    SET_EVAL_B_PHIxPHI(g_flags_eval_id);
    // integrated dipoles 
    USET_FLAG(F_EVAL_D_GG_ALL,g_flags_eval_id);
    SET_FLAG(F_EVAL_D_QG_CONT,g_flags_eval_id);
    
    cout << endl << " EVAL_V_FLAGS = " << bitset<32>(g_flags_eval_v ).to_string() << endl;
    cout << endl << " EVAL_D_FLAGS = " << bitset<32>(g_flags_eval_id).to_string() << endl;

    
    double res_b=0,res_v=0,res_d=0,res_dx=0,res_dx_gg=0,res_dx_qg=0;
    double y_values[] = {-0.9};//-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9};
    double x_values[] = {0.85};//,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};

    for (auto x : x_values)
      {
	for (auto y : y_values)
	  {
	    if (ps_gg_tt.set(sqrt(S),y) )
	      {
	        // set also boosted phase space for integrated dipoles (cont. part)
		ps_gg_tt_x.set(sqrt(x*S),y);
		
		FV const& p1 = ps_gg_tt.p1();
		FV const& p2 = ps_gg_tt.p2();
		FV const& k1 = ps_gg_tt.k1();
		FV const& k2 = ps_gg_tt.k2();
		
#ifdef WITH_T_SPIN
		cout << endl << "Using spin dependent matrix elements. Spin 4-vectors are:" << endl;
		// ALWAYS use spin vectors with the properties
		// s_i.s_i = -1 and s_i.k_i = 0 !!!
		// otherwise tests will fail, since these eqs. are used in the MEs
		FV& s1 = ps_gg_tt.s1();
		FV& s2 = ps_gg_tt.s2();
		//////////////////////////////////////
		// spin vectors in t/tbar restframe
		s1 = {0.0,0.0,sqrt(0.5),+sqrt(0.5)};
		s2 = {0.0,sqrt(0.5),0.0,-sqrt(0.5)};
		// boost from k1 restframe to tt z.m.f.
		static LT boost_k1_RF;
		boost_k1_RF.set_boost(k1,1);
		// boost from k2 restframe to tt z.m.f.
		static LT boost_k2_RF;
		boost_k2_RF.set_boost(k2,1);
		boost_k1_RF.apply(s1);
		boost_k2_RF.apply(s2);
		PRINT(sp(k1,s1));
		PRINT(sp(k2,s2));
		PRINT_4VEC(s1);
		PRINT_4VEC(s2);
		//////////////////////////////////////
#endif



		
		cout << endl;
		// cout << endl << " Setting '" << ps_gg_tt.d_name << "' to ";
		cout << endl << "sqrt(s) = " << ps_gg_tt.get_rs() << "*mScale , y = " << y << ", x = " << x;
		cout << endl << "---------------------------------------------------";
		cout << endl << " Born                   = " << 0.25*(res_b = PREF_GG*Eval_B   (ps_gg_tt,THDM,g_flags_eval_v,GGH_EFF)) << "  " << endl;
		cout << endl << " Born (QQ)              = " << 0.25*(        PREF_QQ*Eval_B_QQ(ps_gg_tt,THDM)) << "  " << endl;		
		cout << endl << " Virt.                  = " << 0.25*(res_v = PREF_GG*Eval_V   (ps_gg_tt,THDM,g_flags_eval_v)) << "  " << endl;
		cout << endl << " DIP                    = " << 0.25*(res_d = PREF_GG*Eval_ID_GG  (ps_gg_tt,x,THDM,g_flags_eval_id)) << "  " << endl;
		Eval_ID_X(ps_gg_tt_x,x,THDM,g_flags_eval_id,res_dx_gg,res_dx_qg);
		cout << endl << " DIP_X                  = " << 0.25*(res_dx= PREF_GG*(res_dx_gg+res_dx_qg)) << "  " << endl;
		cout << endl << " DIP Sum                = " << res_d+res_dx;
		cout << endl << "---------------------------------------------------";
		cout << endl << " NLO                    = " << (res_d+res_dx+res_v+res_b) << "  " << endl;
	      }
	  }
	cout << endl << " =========================================================================================" << endl;
      }
  }

  
  //exit(1);
  
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
	g_flags_eval_r = 0;
	SET_EVAL_R_PHIxPHI_ISR(g_flags_eval_r);
	SET_EVAL_R_PHIxPHI_FSR(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<32>(g_flags_eval_r ).to_string() << endl;
	cout << endl << " 2*RE[M_QCD * M_phi ]     = " << 0.25*(res_r_int = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip 2*RE[M_QCD * M_phi ] = " << 0.25*(res_d_int = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	// phi^2 terms [FSR]
	SET_FLAG(F_EVAL_R_GG,g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_ISR(g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_FSR(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<32>(g_flags_eval_r ).to_string() << " [FSR] " << endl;
	cout << endl << " [M_phi]^2            = " << 0.25*(res_r_phi = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << 0.25*(res_d_phi = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	sum_r += res_r_phi;
	sum_d += res_d_phi;
	// phi^2 terms [ISR]
	USET_EVAL_R_ALL(g_flags_eval_r);
	SET_EVAL_R_ISR_ISR(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<32>(g_flags_eval_r ).to_string() << " [ISR] " << endl;
	cout << endl << " [M_phi]^2            = " << 0.25*(res_r_phi = PREF_GG*Eval_R_GG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << 0.25*(res_d_phi = PREF_GG*Eval_UID_GG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	sum_r += res_r_phi;
	sum_d += res_d_phi;

	// qq
	USET_EVAL_R_ALL(g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_QQ(g_flags_eval_r);
	SET_EVAL_R_PHIxQCD_QQ(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<32>(g_flags_eval_r ).to_string() << " [QQ] " << endl;
	cout << endl << " [M_phi]^2            = " << 0.25*(res_r_phi = PREF_QQ*Eval_R_QQ(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;

	// qg
	USET_EVAL_R_ALL(g_flags_eval_r);
	USET_EVAL_R_PHIxPHI_QG(g_flags_eval_r);
	SET_EVAL_R_PHIxQCD_QG(g_flags_eval_r);
	cout << endl << " EVAL_R_FLAGS = " << bitset<32>(g_flags_eval_r ).to_string() << " [QG] " << endl;
	cout << endl << " [M_phi]^2            = " << 0.25*(res_r_phi = PREF_QG*Eval_R_QG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << 0.25*(res_d_phi = PREF_GG*Eval_UID_QG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;

	cout << endl << " === SWAP !!! === " << endl;
        ps_gg_ttg.swap_initial_state();
	cout << endl << " [M_phi]^2            = " << 0.25*(res_r_phi = PREF_QG*Eval_R_QG(ps_gg_ttg,THDM,g_flags_eval_r))/mScale2 << " [Gev^-2]  " << endl;
	cout << endl << " dip [M_phi]^2        = " << 0.25*(res_d_phi = PREF_GG*Eval_UID_QG(ps_gg_ttg,THDM,g_flags_eval_r)) /mScale2 << " [Gev^-2]  " << endl;
	ps_gg_ttg.swap_initial_state();

	cout << endl << " sum_r = " << sum_r/mScale2;
	cout << endl << " sum_d = " << sum_d/mScale2 << endl;
      }
  }
  return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
