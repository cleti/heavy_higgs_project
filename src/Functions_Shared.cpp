


#include "../inc/Functions_Shared.h"

//#define DEBUG_LT
#define DEBUG_PS_2_3


#define EPS_X     1e-10
#define EPS_BETA2 1e-10
#define EPS_D3    1e-10
#define EPS_D4    1e-10



///// ARG /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/* The invocation name of this program.  */
char *program_name;

struct opt g_options = {
  // default values
  63,//int int_flags;
  10000000,//int n_calls;
  0.0,//double ren_scale;
  14.0,// double cme
  500.0,//double mH;
  3.245439053359254e+01,//double GammaH; // this is the value I use to compare stuff with P.G.
  1.0,//double At;
  1.0,//double Bt;
  0.0,//double Ab;
  0.0,//double Bb;
  7,//int tech_cut;
  6,//int precision; // output
  false,//bool dist;
  false,//bool tedcay;
  false,//logfile
  false,//rootfile
  0//verb_level;
};

struct option const longopts[] = {
  {"int_flags", required_argument, NULL, 'I'},
  {"n_calls", required_argument, NULL, 'N'},
  {"ren_scale", required_argument, NULL, 'R'},
  {"cme", required_argument, NULL, 'E'},
  {"add_higgs", required_argument, NULL, 'H'},
  {"tech_cut", required_argument, NULL, 'T'},
  {"precision", required_argument, NULL, 'P'},
  {"dist", no_argument, NULL, 'D'},
  {"t_decay", no_argument, NULL, 't'},
  {"logfile", no_argument, NULL, 'L'},
  {"rootfile", no_argument, NULL, 'F'},
  {"verbose", required_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"version", no_argument, NULL, 'v'},
  {NULL, 0, NULL, 0}
};







void
usage (int status)
{
  if (status != 0)
    fprintf (stderr, ("Try `%s --help' for more information.\n"),
	     program_name);
  else
    {
      printf (("\
Usage: %s [OPTION]... \n\
 "),
	      program_name);
      fputs (("\
Calculate total and differential cross-section pp -> tt + X at NLO QCD.\n\
\n\
"), stdout);
      fputs (("\
Optional arguments: \n\
"), stdout);
      fputs (("\
  --flags_int,      -I, required argument    integration flags [RDVBB], default: 11111\n\
  --n_calls,        -N, required argument    #integrand calls from VEGAS, default: 1e6\n\
  --ren_scale,      -R, required argument    scale in units of mt,  default: 1.0\n\
  --cme             -E, required argument    pp center of mass energy in TeV, default 14\n\
  --add_higgs       -H  required argument    add higgs boson M,G,a_t,b_t,a_b,b_b [no whitespace!]\n\
  --tech_cut,       -T, required argument    cut on collinear/soft phase space in units of 10^-1,  default: 1.0e-8\n\
  --precision,      -P, no argument          #digits in ouput of numbers, default 6\n\
  --plot,           -D, no argument          plot distributions, default: 0\n\
  --t_decay         -t, no argument          include top/antitop decays and lepton distributions\n\
  --logfile         -L, no argument          write output to file in folder results/\n\
  --rootfile        -F, no argument          store histograms in a ROOT-file in folder results/\n\
  --verbose         -V, required argument    verbosity, default: 0\n \
  --help,           -h, no argument\n\
  --version,        -v, no argument\n\
\n"), stdout);
    }
  exit (status);
}

void version(int status)
{
  printf ("\n VERSION: %s compiled %s : %s \n",program_name,__DATE__,__TIME__);
  exit(status);
}



void parse_arguments(int argc, char** argv, struct opt& options,HiggsModel& hm)
{
  int c;
  std::string arg;
  int pos0,pos1;
  
  while ((c = getopt_long (argc, argv, "I:N:R:E:H:T:P:DtLFV:hv", longopts, NULL)) != -1)
    {
      double m=0.0,g=0.0,at=1.0,bt=1.0,ab=0.0,bb=0.0;
	
      switch (c)
	{
	case 'I':
	  options.int_flags  = std::bitset<16>(std::string(optarg)).to_ulong();
	  break;
	case 'N':
	  options.n_calls = atoi(optarg);
	  break;
	case 'R':
	  options.ren_scale = fabs(atof(optarg))*hm.Scale();
	  break;
	case 'E':
	  options.cme = fabs(atof(optarg));
	  break;
	case 'H':
	  arg = optarg;
	  pos0 = 0;
	  pos1 = 0;

	  // first number is the mass
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      m = atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }
	  else
	    {
	      WARNING("mass is not valid");
	      break;
	    }
	  // second the width
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      g = atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }
	  else
	    {
	      WARNING("width is not valid");
	      break;
	    }
	  // scalar top coupling
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      at = atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }
	  // pseudoscalar top coupling
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      bt =  atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }  
	  // scalar bottom coupling
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      ab = atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }
	  // pseudoscalar bottom coupling
	  if (pos1>=0)
	    {
	      pos1 = arg.find(",",pos0);
	      bb =  atof(arg.substr(pos0,pos1-pos0).c_str());
	      pos0=pos1+1;
	    }

	  hm.AddBoson(m,g,at,bt,ab,bb);
	  break; 
	case 'T':
	  options.tech_cut = abs(atoi(optarg));
	  break;
	case 'P':
	  options.precision = abs(atoi(optarg));
	  break;
	case 'D':
	  options.dist  = !options.dist;
	  break;
	case 't':
	  options.tdecay  = !options.tdecay;
	  break;
	case 'L':
	  options.logfile  = !options.logfile;
	  break;
	case 'F':
	  options.rootfile  = !options.rootfile;
	  break;	  
	case 'V':
	  options.verb_level  = atoi(optarg);
	  break;
	case 'h':
	  usage(0);
	  break;
	case 'v':
	  version(0);
	  break;
	default:
	  usage (1);
	}
    }

  if (!hm.NBosons())
    {
      hm.AddBoson(
		  500.0,
		  3.245439053359254e+01,
		  1.0,
		  1.0,
		  0.0,
		  0.0);
    }
  
  std::cout << std::setprecision(g_options.precision);
}
///// ARG /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////







/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////





CanvasPtr MakeCanvas(const char* title,
		     int width,
		     int height,
		     TCanvas* c,
		     double left_marg_scale,
		     double bottom_marg_scale)
{
  static int c_counter = 0;
  static int p_counter = 0;
  // alllocate new canvas if neccessary
  if ( c == nullptr )
    {
      c = new TCanvas(cat("Canvas_",c_counter++),title,0,0,width,height);
      c->Divide(1,1);
    }
  double can_size_ratio = (double)height/(double)width;

  // upper pad
  TPad* p1_1 = new TPad(cat("upper_pad_",p_counter),cat("upper_pad_",p_counter),0,0.2,1,1);
  p1_1->SetTopMargin(0.05);
  p1_1->SetBottomMargin(0.0075);
  p1_1->SetLeftMargin(0.1*can_size_ratio*left_marg_scale);
  p1_1->SetRightMargin(0.05);
  // lower pad
  TPad* p1_2 = new TPad(cat("lower_pad_",p_counter),cat("lower_pad_",p_counter),0,0,1,0.2);
  p1_2->SetTopMargin(0.0025);
  p1_2->SetBottomMargin(0.3*bottom_marg_scale);
  p1_2->SetLeftMargin(0.1*can_size_ratio*left_marg_scale);
  p1_2->SetRightMargin(0.05);
  return {c,p1_1,p1_2};
}
DoubleCanvasPtr MakeDoubleCanvas(const char* title, int width, int height)
{
  static int c_counter = 0;
  TCanvas* c = new TCanvas(cat("DoubleCanvas_",c_counter++),title,0,0,width,height);
  c->Divide(2,1);
  return {MakeCanvas(title,width/2,height,c),MakeCanvas(title,width/2,height,c)};
}


// sum weights in the histogram
static double GetIntegral(TH1D* hist)
{
  double sum = 0.0;
  int n = hist->GetNbinsX();
  for (int i=1;i<=n;++i)
    {
      sum += (hist->GetBinContent(i))*(hist->GetBinWidth(i));
    }
  return sum;
}


void DrawDistribution(CanvasPtr& canvas,
		      HistArray& histograms,
		      double* norm,
		      std::string const& titleX,
		      std::string const& titleY,
		      bool WRITE,
		      bool NLO,
		      int CD)
{
  double TitleOffsetY = 1.2;
  double LabelSizeY   = 0.028;
  double LegXPos      = 0.6557;
  // single plot
  if (CD==0)
    {
      TitleOffsetY = 0.8;
      LegXPos = 0.7428;
    }
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  TLegend* leg = new TLegend(
			     LegXPos,
			     0.8461 ,
			     0.9559 ,
			     0.9550
			     );

      canvas.c->cd(CD); // upper pad conatins the distributions
      canvas.p1_1->Draw();
      canvas.p1_1->cd();

      // 0 contains LO QCD
      histograms[0]->Smooth(1);
      // 1 conatins LO QCD + LO PHI^2 + LO PHIxQCD
      histograms[1]->Add(histograms[0]);
      if (NLO)
	{
	  // 3 conatins complete NLO
	  histograms[3]->Add(histograms[1]);
	  // 4 conatins NLO virt.
	  histograms[3]->Add(histograms[4]);
	  // 5 contains NLO real
	  histograms[5]->Smooth(1);
	  histograms[3]->Add(histograms[5]);
	  
	}

      // normalize
      histograms[0]->Scale(1.0/norm[0]);// LO QCD cross section
      histograms[1]->Scale(1.0/norm[1]);// LO QCD+PHI cross section

      if (NLO)
	{
	  histograms[3]->Scale(1.0/norm[3]);// NLO QCD+PHI cross section

	  histograms[3]->SetLineWidth(2);
	  histograms[3]->SetLineColor(2);
	  histograms[3]->GetYaxis()->SetLabelSize(LabelSizeY);
	  histograms[3]->GetYaxis()->SetTitleOffset(TitleOffsetY);
	  histograms[3]->GetYaxis()->SetTitle(titleY.c_str());
	  histograms[3]->DrawCopy("hist ][");

	  histograms[1]->SetLineWidth(1);
	  histograms[1]->SetLineColor(2);
	  histograms[1]->DrawCopy("same hist ][");
	}
      else
	{
	  histograms[1]->SetLineWidth(2);
	  histograms[1]->SetLineColor(2);
	  histograms[1]->GetYaxis()->SetLabelSize(LabelSizeY);
	  histograms[1]->GetYaxis()->SetTitleOffset(TitleOffsetY);
	  histograms[1]->GetYaxis()->SetTitle(titleY.c_str());
	  histograms[1]->DrawCopy("hist ][");
	}

      histograms[0]->SetLineWidth(1);
      histograms[0]->SetLineColor(1);
      histograms[0]->DrawCopy("same hist ][");
      
      leg->AddEntry(histograms[0] ,"LO QCD","L");
      if (RunParameters::TwoHDM)
	{
	  leg->AddEntry(histograms[1] ,"LO QCD + #Phi_{1} + #Phi_{2} ","L");
	}
      else
	{
	  leg->AddEntry(histograms[1] ,"LO QCD + #Phi ","L");
	}
      if (NLO)
	{
	  if (RunParameters::TwoHDM)
	    {
	      leg->AddEntry(histograms[3] ,"NLO QCD + #Phi_{1} + #Phi_{2} ","L");
	    }
	  else
	    {
	      leg->AddEntry(histograms[3] ,"NLO QCD + #Phi ","L");
	    }
	}
      leg->Draw("same");

      //////////////////////////////////////////////////////
      std::cout << std::endl;
      std::cout << histograms[0]->GetName() << std::endl;
      std::cout << " I[LO]   = " << GetIntegral(histograms[0]) << std::endl;
      std::cout << " I[LO*]  = " << GetIntegral(histograms[1]) << std::endl;
      std::cout << " I[NLO]  = " << GetIntegral(histograms[3]) << std::endl;
      //////////////////////////////////////////////////////
      canvas.c->cd(CD); //lower pad contains ratio NLO/LO
      canvas.p1_2->Draw();
      canvas.p1_2->cd();

      histograms[3]->Divide(histograms[0]);
      histograms[1]->Divide(histograms[0]);
      histograms[0]->Divide(histograms[0]);

      if (NLO)
	{
	  SetRatioPlot(histograms[3],titleX.c_str());
	  histograms[3]->DrawCopy("hist ][");
	  histograms[1]->DrawCopy("same hist ][");
	}
      else
	{
	  SetRatioPlot(histograms[1],titleX.c_str());
	  histograms[1]->DrawCopy("hist ][");
	}
      histograms[0]->DrawCopy("same hist ][");
      //////////////////////////////////////////////////////
      //////////////////////////////////////////////////////
      if (WRITE) canvas.c->Write();
}


void DrawTwoDistributions(DoubleCanvasPtr& canvas,
			  HistArray& histogramsL,
			  HistArray& histogramsR,
			  double* norm,
			  std::string const& titleLX,std::string const& titleRX,
			  std::string const& titleLY,std::string const& titleRY,
			  bool WRITE,
			  bool NLO)
{
  DrawDistribution(canvas.c1,
		   histogramsL,
		   norm,
		   titleLX,
		   titleLY,
		   WRITE,
		   NLO,
		   1);
  DrawDistribution(canvas.c2,
		   histogramsR,
		   norm,
		   titleRX,
		   titleRY,
		   WRITE,
		   NLO,
		   2);
}

void SetRatioPlot(TH1D* hist,std::string xtitle)
{
  hist->SetTitle("");
  hist->Smooth();
  hist->GetYaxis()->SetTitleSize(12); //in pixels
  hist->GetYaxis()->SetLabelFont(43); //font in pixels
  hist->GetYaxis()->SetLabelSize(18); //in pixels
  hist->GetYaxis()->SetTitleOffset(3.5);
  hist->GetYaxis()->SetTitle("NLO / LO");
  hist->GetYaxis()->SetTickLength(0.02); //in pixels
  hist->GetYaxis()->SetNdivisions(505,1); //in pixels
  //hist->GetYaxis()->SetRangeUser(-0.29,0.29);
  hist->GetYaxis()->SetRangeUser(0.91,1.09);
  hist->GetXaxis()->SetTitleFont(43); //font in pixels
  hist->GetXaxis()->SetTitleSize(12); //in pixels
  hist->GetXaxis()->SetTitleOffset(3.5);
  hist->GetXaxis()->SetLabelFont(43); //font in pixels
  hist->GetXaxis()->SetLabelSize(20); //in pixels
  hist->GetXaxis()->SetTickLength(0.1); //in pixels
  hist->GetXaxis()->SetTitle(xtitle.c_str());
  hist->GetXaxis()->SetTitleSize(25);
  hist->GetXaxis()->SetTitleOffset(4.5);
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////































////////////////////////////////////////////////////////////////////////////
// complex Phi propagator denominator
////////////////////////////////////////////////////////////////////////////
namespace NSP_DENS
{
  double   den_t  = 0.0;
  c_double cden_t = c_double(0.0,0.0);
}
// abs. value squared
double const& DenS2(double const& s, double const& m, double const& g)
{
  using namespace NSP_DENS;
  double const& m2 = m*m;
  double const& mg = m*g;
  den_t  = 0.1e1 / (pow(s-m2,2)+pow(mg,2));
  // important!
  // if this is not here, DenS will return 0 if DenS2 gets called first and s is unchanged!
  cden_t = c_double((s-m2)*den_t,-mg*den_t);
  return den_t;
}
// real and imaginary parts in a complex double
c_double const& DenS(double const& s, double const& m, double const& g)
{
  using namespace NSP_DENS;
  double const& m2 = m*m;
  double const& mg = m*g;
  den_t  = 0.1e1 / (pow(s-m2,2)+pow(mg,2));
  cden_t = c_double((s-m2)*den_t,-mg*den_t);
  return cden_t;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


// this function computes the spin 4-vectors s1,s2 in the k1,k2-z.m.f. (given k1 and k2 in precisely this frame)
// assuming that s1_r and s2_r are the spin directions in the respective restframes
// the restframes are connected to the k1,k2-z.m.f. by rotationfree boosts in the respective directions
// the spin 4-vectors then obey k_i.s_i = 0 and s_i.s_i = -1
void set_spins_in_tt_zmf(
			 FV const& k1,
			 FV const& k2,
			 FV& s1,
			 FV& s2,
			 FV const& s1_r,
			 FV const& s2_r
			 )
{
  double e1   = (1.0+k1[0]);
  double k1s1 = (k1[1]*s1_r[1]+k1[2]*s1_r[2]+k1[3]*s1_r[3])/e1;
  s1[0] = k1s1*e1;
  s1[1] = s1_r[1]+k1[1]*k1s1;
  s1[2] = s1_r[2]+k1[2]*k1s1;
  s1[3] = s1_r[3]+k1[3]*k1s1;
  double e2   = (1.0+k2[0]);
  double k2s2 = (k2[1]*s2_r[1]+k2[2]*s2_r[2]+k2[3]*s2_r[3])/e2;
  s2[0] = k2s2*e2;
  s2[1] = s2_r[1]+k2[1]*k2s2;
  s2[2] = s2_r[2]+k2[2]*k2s2;
  s2[3] = s2_r[3]+k2[3]*k2s2;
}





struct eps_entry {
  int  indices[4];
  int  sign;
};

// double EPS_(double const* __restrict__ k1, double const* __restrict__ k2, double const* __restrict__ k3, double const* __restrict__ k4, int DEALLOC)
double EPS_(FV const& k1, FV const& k2, FV const& k3, FV const& k4)
{
  /////////////////////////////////////////////////////////////////////////////////////////////
  // executed only on first call //////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  static bool DO_INIT = 1;
  // 4dim epsilon tensor entries in linear row
  static eps_entry eps_list[24];
  /////////////////////////////////////////////////////////////////////////////////////////////
  if ((DO_INIT)) {
    // initialize eps_list entries
    double eps[4][4][4][4];
    int I = 0;
    for (int i=0;i<4;++i)
      {
	for (int j=0;j<4;++j)
	  {
	    for (int k=0;k<4;++k)
	      {
		for (int l=0;l<4;++l)
		  {
		    eps[i][j][k][l] = (double)((i-j)*(i-k)*(i-l)*(j-k)*(j-l)*(k-l)/12);
		    if (eps[i][j][k][l]!=0.0)
		      {
			eps_list[I].indices[0] = i;
			eps_list[I].indices[1] = j;
			eps_list[I].indices[2] = k;
			eps_list[I].indices[3] = l;
			eps_list[I++].sign = eps[i][j][k][l];
		      }
		  }
	      }
	  }
      }
    DO_INIT = 0;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  // this code is executed each call
  // loop is reduced to the 24 non-zero elements in the eps-tensor
  /////////////////////////////////////////////////////////////////////////////////////////////
  double res = 0.0;
  for (int I=0; I<24; ++I)
    {
      int const& i    = eps_list[I].indices[0];
      int const& j    = eps_list[I].indices[1];
      int const& k    = eps_list[I].indices[2];
      int const& l    = eps_list[I].indices[3];
      if (k1[i]==0.0 || k2[j]==0.0 || k3[k]==0.0 || k4[l]==0.0) continue;
      res -= k1[i]*k2[j]*k3[k]*k4[l]*eps_list[I].sign;
    }
  return res;
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
}



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
namespace AmpPrefactors {

  void ResetAmpPrefactors()
  {
    using namespace Constants;
    using namespace RunParameters;
    /* born */
    PREF_B_PHIxPHI      = CF*CA2*AlphaS2/Pi2;
    PREF_B_PHIxQCD      = CF*CA*AlphaS2;
    PREF_B_QCDxQCD      = CF*CA*AlphaS2*Pi2;
    PREF_B_QCDxQCD_CF   = CF2*CA*AlphaS2*Pi2;
    PREF_B_QCDxQCD_CA   = CF*CA2*AlphaS2*Pi2;
    PREF_B_QCDxQCD_CFCA = CF*CA*CFCA2*AlphaS2*Pi2;
    /* virtual */
    PREF_V_CF     = PREF_B_PHIxQCD*AlphaS/Pi*CF;
    PREF_V_CA     = PREF_B_PHIxQCD*AlphaS/Pi*CA;
    PREF_V_CFCA2  = PREF_B_PHIxQCD*AlphaS/Pi*CFCA2;
    PREF_V_Nf     = PREF_B_PHIxQCD*AlphaS/Pi*Nf;
    PREF_V_CT     = PREF_B_PHIxQCD*AlphaS/Pi;
    PREF_V        = PREF_B_PHIxQCD*AlphaS;
    PREF_V_PHI    = PREF_B_PHIxPHI*AlphaS/Pi; // CF*CA2*AlphaS^3/Pi^3
    PREF_V_PHI_CA = PREF_B_PHIxPHI*AlphaS/Pi*CA;
    PREF_V_PHI_CF = PREF_B_PHIxPHI*AlphaS/Pi*CF;
    /* real */
    PREF_R        = PREF_V_CT*Pi2;
    PREF_R_CF     = PREF_V_CF*Pi2;
    PREF_R_CA     = PREF_V_CA*Pi2;
    PREF_R_CFCA2  = PREF_V_CFCA2*Pi2;
    PREF_R_PHI    = PREF_V_PHI*Pi2; // CF*CA2*AlphaS^3/Pi
    PREF_R_PHI_CA = PREF_V_PHI_CA*Pi2;
    PREF_R_PHI_CF = PREF_V_PHI_CF*Pi2;
    /* UID */
    PREF_UID_TF =  4.0*Pi*AlphaS;
    PREF_UID_CA = 16.0*Pi*AlphaS*CA;
    PREF_UID_CF =  8.0*Pi*AlphaS*CF;
  }
  void PrintAmpPrefactors()
  {
    std::cout << "\n\n Amplitude-Prefactors:\n";
    /* born */
    PRINT(    PREF_B_PHIxPHI     );
    PRINT(    PREF_B_PHIxQCD     );
    PRINT(    PREF_B_QCDxQCD     );
    PRINT(    PREF_B_QCDxQCD_CF  );
    PRINT(    PREF_B_QCDxQCD_CA  );
    PRINT(    PREF_B_QCDxQCD_CFCA);
    /* virtual */
    PRINT(    PREF_V_CF     );
    PRINT(    PREF_V_CA     );
    PRINT(    PREF_V_CFCA2  );
    PRINT(    PREF_V_Nf     );
    PRINT(    PREF_V_CT     );
    PRINT(    PREF_V        );
    PRINT(    PREF_V_PHI_CA );
    PRINT(    PREF_V_PHI_CF );
    /* real */
    PRINT(    PREF_R       );
    PRINT(    PREF_R_CF    );
    PRINT(    PREF_R_CA    );
    PRINT(    PREF_R_CFCA2 );
    PRINT(    PREF_R_PHI_CA);
    PRINT(    PREF_R_PHI_CF);
    /* UID */
    PRINT(    PREF_UID_CA);
    PRINT(    PREF_UID_CF);
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace HiggsBosons {
  void SetMPhi1(double const& M)
  {
    using namespace RunParameters;
    M_1  = M/mScale;
    M2_1 = M_1*M_1;
  }
  // setter function Phi1 (M,Gamma in [GeV])
  void SetPhi1(double const& M,
	       double const& Gamma,
	       double const& at,
	       double const& bt,
	       double const& ab,
	       double const& bb)
  {
    using namespace RunParameters;
    M_1  = M/mScale;
    M2_1 = M_1*M_1;
    G_1  = Gamma/mScale;
    G2_1 = G_1*G_1;
    At_1 = at/Vh;
    Bt_1 = bt/Vh;
    Ab_1 = ab/Vh;
    Bb_1 = bb/Vh;
    FH_eff_1 = (At_1)/12.0;
    FA_eff_1 = (bt!=0)?(-Bt_1/16.0):0.0;
  }
  // setter function Phi2 (M,Gamma in [GeV])
  void SetPhi2(double const& M,
	       double const& Gamma,
	       double const& at,
	       double const& bt,
	       double const& ab,
	       double const& bb)
  {
    using namespace RunParameters;
    M_2  = M/mScale;
    M2_2 = M_2*M_2;
    G_2  = Gamma/mScale;
    G2_2 = G_2*G_2;
    At_2 = at/Vh;
    Bt_2 = bt/Vh;
    Ab_2 = ab/Vh;
    Bb_2 = bb/Vh;
    FH_eff_2 =  At_2/12.0;
    FA_eff_2 = (bt!=0)?(-Bt_2/16.0):0.0;
  }
  void PrintHiggsPrefactors()
  {
    std::cout << "\n\n Higgsboson-specific prefactors:\n";
    // interference terms
    PRINT(At_fH_re);		
    PRINT(At_fA_re);		
    PRINT(Bt_fH_re);		
    PRINT(Bt_fA_re);
  
    PRINT(At_fH_im);
    PRINT(At_fA_im);	
    PRINT(Bt_fH_im);		
    PRINT(Bt_fA_im);
    // phi^2
    PRINT(At2_fH2_De);	
    PRINT(At2_fA2_De);
    PRINT(Bt2_fH2_De);	
    PRINT(Bt2_fA2_De);
#ifdef WITH_T_SPIN
    PRINT(At_Bt_fH2_De);	
    PRINT(At_Bt_fA2_De);
#endif
  }
  
  void ResetHiggsPrefactors(double const& S, unsigned EFF)
  {
    // static double S_t = 0.0;
    // if (S != S_t)
    //   {
    // 	S_t = S;
	c_double FH_1 = 0.0;
	c_double FA_1 = 0.0;
	
	switch (EFF)
	  {
	  case 0: // full ggH vertex
	    set_F_ggH(S);
	    FH_1 = -F_ggH1_s/(4.0*S);
	    FA_1 =  F_ggH1_p/(8.0*S);
	    break;
	  case 1: // effective ggH vertex
	    FH_1 = FH_eff_1;
	    FA_1 = FA_eff_1;
	    break;
	  }
	
	// PRINT(RunParameters::AlphaS/Constants::Pi*FH_1);
	// PRINT(RunParameters::AlphaS/Constants::Pi*FA_1);
	
	double FH2_1  = std::norm(FH_1);
	double FA2_1  = std::norm(FA_1);

	// Phi1 propagator denominator
	double   d_1 = 1.0 / (pow(S-M2_1,2)+pow(M_1*G_1,2));
	c_double D_1 = c_double((S-M2_1)*d_1,-M_1*G_1*d_1);
#ifdef DEBUG
	CHECKNA(d_1);
	CHECKNA(D_1.real());
	CHECKNA(D_1.imag());
#endif

	// PRINT(D_1);
	// PRINT(D_1);

	// PRINT(RunParameters::AlphaS/Constants::Pi*FH_1*D_1);
	// PRINT(RunParameters::AlphaS/Constants::Pi*FA_1*D_1);
	
	// PHIxQCD amplitudes: add Phi1 contribution
	At_fH_re = At_1*(FH_1*D_1).real();
	At_fA_re = At_1*(FA_1*D_1).real();
	Bt_fH_re = Bt_1*(FH_1*D_1).real();
	Bt_fA_re = Bt_1*(FA_1*D_1).real();

	At_fH_im = At_1*(FH_1*D_1).imag();
	At_fA_im = At_1*(FA_1*D_1).imag();
	Bt_fH_im = Bt_1*(FH_1*D_1).imag();
	Bt_fA_im = Bt_1*(FA_1*D_1).imag();

	double At2_1 = At_1*At_1;
	double Bt2_1 = Bt_1*Bt_1;
	
	// PHIxPHI amplitudes: add Phi1^2 contribution
	At2_fH2_De   = At2_1*FH2_1*d_1;
	At2_fA2_De   = At2_1*FA2_1*d_1;
	Bt2_fH2_De   = Bt2_1*FH2_1*d_1;
	Bt2_fA2_De   = Bt2_1*FA2_1*d_1;
	
#ifdef WITH_T_SPIN
	// only needed for polarized amplitudes
	// PHIxPHI terms: add at*bt interference from Phi1^2 contribution
	At_Bt_fH2_De = At_1*Bt_1*FH2_1*d_1;
	At_Bt_fA2_De = At_1*Bt_1*FA2_1*d_1;
#endif
	
	// add the contribution of the second boson if required
	if (RunParameters::TwoHDM>0)
	  {
	    c_double FH_2 = 0.0;
	    c_double FA_2 = 0.0;
	
	    switch (EFF)
	      {
	      case 0: // full ggH vertex
		FH_2 = -F_ggH2_s/(4.0*S);
		FA_2 =  F_ggH2_p/(8.0*S);
		break;
	      case 1: // effective ggH vertex
		FH_2 = FH_eff_2;
		FA_2 = FA_eff_2;
		break;
	      }

	    double FH2_2  = std::norm(FH_2);
	    double FA2_2  = std::norm(FA_2);
	    
	    // Phi2 propagator denominator
	    double   d_2 = 1.0 / (pow(S-M2_2,2)+pow(M_2*G_2,2));
	    c_double D_2 = c_double((S-M2_2)*d_2,-M_2*G_2*d_2);
#ifdef DEBUG
	    CHECKNA(d_2);
	    CHECKNA(D_2.real());
	    CHECKNA(D_2.imag());
#endif

	    // PHIxQCD amplitudes: add Phi2 contribution
	    At_fH_re += At_2*(FH_2*D_2).real();
	    At_fA_re += At_2*(FA_2*D_2).real();
	    Bt_fH_re += Bt_2*(FH_2*D_2).real();
	    Bt_fA_re += Bt_2*(FA_2*D_2).real();

	    At_fH_im += At_2*(FH_2*D_2).imag();
	    At_fA_im += At_2*(FA_2*D_2).imag();
	    Bt_fH_im += Bt_2*(FH_2*D_2).imag();
	    Bt_fA_im += Bt_2*(FA_2*D_2).imag();
	    
	    double At2_2 = At_2*At_2;
	    double Bt2_2 = Bt_2*Bt_2;
	    
	    // PHIxPHI amplitudes: add Phi2^2 contribution
	    At2_fH2_De   += At2_2*FH2_2*d_2;
	    At2_fA2_De   += At2_2*FA2_2*d_2;
	    Bt2_fH2_De   += Bt2_2*FH2_2*d_2;
	    Bt2_fA2_De   += Bt2_2*FA2_2*d_2;

#ifdef WITH_T_SPIN
	    // only needed for polarized amplitudes
	    // PHIxPHI terms: add at*bt interference from Phi1^2 contribution
	    At_Bt_fH2_De += At_2*Bt_2*FH2_2*d_2;
	    At_Bt_fA2_De += At_2*Bt_2*FA2_2*d_2;
#endif
	    
	    if (RunParameters::TwoHDM==1)
	      {
		// PHIxPHI terms: Phi1 * Phi2 interferences
		c_double D1D2 = D_1*std::conj(D_2);
		// boson factors for phi^2: Phi1-Phi2 interference terms
		At2_fH2_De   += 2.0*At_1*At_2*(FH_1*std::conj(FH_2)*D1D2).real();
		At2_fA2_De   += 2.0*At_1*At_2*(FA_1*std::conj(FA_2)*D1D2).real();
		Bt2_fH2_De   += 2.0*Bt_1*Bt_2*(FH_1*std::conj(FH_2)*D1D2).real();
		Bt2_fA2_De   += 2.0*Bt_1*Bt_2*(FA_1*std::conj(FA_2)*D1D2).real();
#ifdef WITH_T_SPIN
		// only needed for polarized amplitudes
		// PHIxPHI terms: at Phi1 * bt  Phi2 interferences
		At_Bt_fH2_De += (At_1*Bt_2+At_2*Bt_1)*(FH_1*FH_2*D1D2).real();
		At_Bt_fA2_De += (At_1*Bt_2+At_2*Bt_1)*(FA_1*FA_2*D1D2).real();	
		// the Phi1 At * Phi2 Bt interference is prop. to the imaginary parts
		At_Bt_fH2_DeIM = (At_1*Bt_2-At_2*Bt_1)*(FH_1*std::conj(FH_2)*D1D2).imag();
		At_Bt_fA2_DeIM = (At_1*Bt_2-At_2*Bt_1)*(FA_1*std::conj(FA_2)*D1D2).imag();
#endif
	      }	
	  }
      // }
  }


 // static c_double I3_test(double S, double M2)
 //  {
 //    double tau = 4.0*M2/S;
 //    if (tau>=1.0)
 //      {
 // 	return -2.0/S*pow(std::asin(1.0/sqrt(tau)),2);
 //      }
 //    else
 //      {
 // 	double beta = sqrt(1.0-tau);
 // 	return 1.0/(2.0*S)*std::pow(c_double(std::log((1.0+beta)/(1.0-beta)),-Constants::Pi),2);
 //      }
 //  }
  
  // set scalar and pseudoscalar ggH form factors (LO)
  // b-quark effects included if the couplings Ab_1,Bb_1 are set
  void set_F_ggH(double S)
  {
    using namespace RunParameters;
    static int I = 0;// finite part (integral here has no divergent part)
    static double zero = 0.0;

    c_double I3t = 0.0;
    c_double I3b = 0.0;
    F_ggH1_s = 0.0;
    F_ggH1_p = 0.0;
	
    // top contribution
    if (At_1 != 0.0 || Bt_1 != 0.0)
      {
	I3t = qli3_(&S,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I);
	// PRINT(I3t);
	// PRINT(I3_test(S,mt2));
	F_ggH1_s += mt2 * At_1 * ( (S-4.0*mt2) * I3t - 2.0 );
	F_ggH1_p += mt2 * Bt_1 * S * I3t;
      }
    // bottom contribution
    if (Ab_1 != 0.0 || Bb_1 != 0.0)
      {
	I3b = qli3_(&S,&zero,&zero,&mb2,&mb2,&mb2,&MUR2,&I);
	// PRINT(I3b);
	// PRINT(I3_test(S,mb2));
	F_ggH1_s += mb2 * Ab_1 * ( (S-4.0*mb2) * I3b - 2.0 );
	F_ggH1_p += mb2 * Bb_1 * S * I3b;
      }
  
    if (TwoHDM>0)
      {
	F_ggH2_s = 0.0;
	F_ggH2_p = 0.0;
	// top contribution
	if (At_2 != 0.0 || Bt_2 != 0.0)
	  {
	    F_ggH2_s += mt2 * At_2 * ( (S-4.0*mt2) * I3t - 2.0 );
	    F_ggH2_p += mt2 * Bt_2 * S * I3t;
	  }
	// bottom contribution
	if (Ab_2 != 0.0 || Bb_2 != 0.0)
	  {
	    F_ggH2_s += mb2 * Ab_2 * ( (S-4.0*mb2) * I3b - 2.0 );
	    F_ggH2_p += mb2 * Bb_2 * S * I3b;
	  }
      }
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace RunParameters {
  void SetAlphaS(double const& val)
  {
    AlphaS  = val;
    AlphaS2 = val*val;
    // needs to be done whenever AlphaS changes
    AmpPrefactors::ResetAmpPrefactors();
  }
  void SetMUR(double const& val)
  {
#ifdef DEBUG
    if (val>10.0) WARNING("this function expects argument in units of mScale");
#endif
    MUR   = val;
    MUR2  = val*val;
    LNMU2 = log(MUR2/mt2);
  }
  void SetMUF(double const& val)
  {
#ifdef DEBUG
    if (val>10.0) WARNING("this function expects argument in units of mScale");
#endif
    MUF  = val;
    MUF2 = val*val;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
