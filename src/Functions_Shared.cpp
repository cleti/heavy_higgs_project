


#include "../inc/Functions_Shared.h"

//#define DEBUG_LT
#define DEBUG_PS_2_3


#define EPS_X     1e-10
#define EPS_BETA2 1e-10
#define EPS_D3    1e-10
#define EPS_D4    1e-10



////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/* The invocation name of this program.  */
char *program_name;

struct opt g_options = {
  // default values
  63,//int int_flags;
  10000000,//int n_calls;
  0.0,//double ren_scale;
  14.0,// double cme
  173.34,//double mt;
  4.75,//double mb; 
  246.0,//double vh;
  7,//int tech_cut;
  6,//int precision; // output
  false,//bool dist;
  false,//bool tdecay;
  false,//bool useK;
  false,//bool logfile;
  false,//bool rootfile;
  0//intverb_level;
};



struct option const longopts[] = {
  {"int_flags", required_argument, NULL, 'I'},
  {"n_calls", required_argument, NULL, 'N'},
  {"ren_scale", required_argument, NULL, 'R'},
  {"cme", required_argument, NULL, 'E'},
  {"mt", required_argument, NULL, 'M'},
  {"mb", required_argument, NULL, 'm'},  
  {"Vh", required_argument, NULL, 'V'},    
  {"add_higgs", required_argument, NULL, 'H'},
  {"tech_cut", required_argument, NULL, 'T'},
  {"precision", required_argument, NULL, 'P'},
  {"verbose", required_argument, NULL, 'v'},
  {"dist", no_argument, NULL, 'D'},
  {"t_decay", no_argument, NULL, 't'},
  {"use_Kfactor", no_argument, NULL, 'K'},
  {"logfile", no_argument, NULL, 'L'},
  {"rootfile", no_argument, NULL, 'F'},
  {"info", no_argument, NULL, 'i'},
  {NULL, 0, NULL, 0}
};







void
usage (int status)
{
    
  if (status != 0)
    fprintf (stderr, ("Try `%s --info' for more information.\n"),
	     program_name);
  else
    {
      printf ("\n Calculate total and differential cross-section pp -> tt + X at NLO QCD.\n");  
      printf ("\n VERSION: %s compiled %s : %s \n",program_name,__DATE__,__TIME__);
      printf ("\n Usage: %s [OPTION]... \n",program_name);

      fputs (("Optional arguments: \n"), stdout);
      fputs (("\
  --flags_int,      -I, required argument    integration flags [RDVBB], default: 11111\n\
  --n_calls,        -N, required argument    #integrand calls from VEGAS, default: 1e6\n\
  --ren_scale,      -R, required argument    scale in units of GeV,  default: 1.0\n\
  --cme             -E, required argument    pp center of mass energy in TeV, default 14\n\
  --mt              -M, required argument    top-quark mass in GeV, default 173.34\n\
  --Vh              -V, required argument    combined Higgs VEV in GeV, default 246\n\
  --add_higgs       -H  required argument    add higgs boson M,G,a_t,b_t,a_b,b_b [no whitespace!]\n\
  --tech_cut,       -T, required argument    cut on collinear/soft phase space in units of 10^-1,  default: 1.0e-8\n\
  --precision,      -P, no argument          #digits in ouput of numbers, default 6\n\
  --plot,           -D, no argument          plot distributions, default: 0\n\
  --t_decay         -t, no argument          include top/antitop decays and lepton distributions\n\
  --use_Kfactor     -K, no argument          rescale effective Higgs-gluon vertex in NLO corrections\n\
  --logfile         -L, no argument          write output to file in folder results/\n\
  --rootfile        -F, no argument          store histograms in a ROOT-file in folder results/\n\
  --verbose         -v, required argument    verbosity level, default: 0\n\
  --info,           -i, no argument\n\
\n"), stdout);
    }
  exit (status);
}



void parse_arguments(int argc, char** argv, struct opt& options,HiggsModel& hm)
{
  const char* arglist = "M:I:N:R:E:m:V:H:T:P:v:DtLFKi";
  std::string arg;

  int pos0=0,pos1=0,c1=0;
  // first scan for all options except Higgs boson parameters
  while ((c1 = getopt_long (argc, argv, arglist, longopts, NULL)) != -1)
    {	
      switch (c1)
	{
	case 'H':
	  break;
	case 'M':
  	  options.mt = fabs(atof(optarg));
  	  break;
	case 'V':
	  options.vh = fabs(atof(optarg));
  	  break;	  
	case 'I':
	  options.int_flags  = std::bitset<16>(std::string(optarg)).to_ulong();
	  break;
	case 'N':
	  options.n_calls = atoi(optarg);
	  break;
	case 'R':
	  options.ren_scale = fabs(atof(optarg));
	  break;
	case 'E':
	  options.cme = fabs(atof(optarg));
	  break;
	case 'm':
	  options.mb = fabs(atof(optarg));
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
	case 'K':
	  options.useK  = !options.useK;
	  break;	  
	case 'L':
	  options.logfile  = !options.logfile;
	  break;
	case 'F':
	  options.rootfile  = !options.rootfile;
	  break;	  
	case 'v':
	  options.verb_level  = atoi(optarg);
	  break;
	case 'i':
	  usage(0);	  
	  break;
	default:
	  usage(1);
	}
    }

  // FIRST set the scale in the HiggsModel
  // from now on all quantities will be normalized to mt
  hm.SetScale(options.mt);
  // set the Higgs VEV
  hm.SetVH(options.vh);
  // set the top-quark mass
  hm.SetMt(options.mt);
  // set the bottom-quark mass
  hm.SetMb(options.mb);
  // use the scaling factor for NLO corrections 
  hm.SetUseK(options.useK);


  // reset the index, so that getopt_long() searches through all the options again
  // (it points to the -M option in argv at this point if the option was provided by the user)
  optind = 0;
  int c2=0;
  // the order is crucial for the extraction of command line arguments:
  // the Higgs bosons can not be set before the scale and VEV are specified
  while ((c2 = getopt_long (argc, argv, arglist, longopts, NULL)) != -1)
    {
      double m=0.0,g=0.0,at=1.0,bt=1.0,ab=0.0,bb=0.0;
      switch (c2)
  	{
	case 'H':
	  arg = optarg;
	  pos0 = 0;
	  pos1 = 0;

	  // scan through the argument string, should have the format
	  // 'M,G,a_t,b_t,a_b,b_b' -- NO WHITESPACE!
	  // first number is the mass
	  if (pos1>=0)
	    {
	      // scan from pos0 and find position of next ',' -> this is pos1
	      pos1 = arg.find(",",pos0);
	      // convert string ranging from pos0 to pos1 into a float
	      m = atof(arg.substr(pos0,pos1-pos0).c_str());
	      // update pos0 -> this is now the first char behind the ','
	      pos0=pos1+1;
	      // and so on ...
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
	}
    }
  


  
  // default boson if no user input via --add_boson ... was provided
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


  
  // set the scale to the mean value of the Higgs masses if no explicit user input was provided
  double MU = options.ren_scale;
  if (MU==0.0)
    {
      int k=0;
      for (auto phi_i = std::begin(hm.GetBosons()); phi_i!=std::end(hm.GetBosons()); ++phi_i)
	{
	  MU += (hm.Scale()*(*phi_i)->M() - MU) / (k+1.0);
	  ++k;
	}
      MU *= 0.5;
    }

  hm.SetMUR(MU);
  hm.SetMUF(MU);
  
  
  std::cout << std::setprecision(options.precision);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////







////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
static std::string cat(std::string str, int number)
{
  std::ostringstream strstr;
  strstr << str << number;
  return strstr.str();
}



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

      c = new TCanvas(cat("Canvas_",c_counter++).c_str(),title,0,0,width,height);
      c->Divide(1,1);
    }
  double can_size_ratio = (double)height/(double)width;

  
  // upper pad
  TPad* p1_1 = new TPad(
			cat("upper_pad_",p_counter).c_str(),
			cat("upper_pad_",p_counter).c_str(),
			0,0.2,1,1);
  ++p_counter;
  p1_1->SetTopMargin(0.05);
  p1_1->SetBottomMargin(0.0075);
  p1_1->SetLeftMargin(0.1*can_size_ratio*left_marg_scale);
  p1_1->SetRightMargin(0.05);
  // lower pad
  TPad* p1_2 = new TPad(
			cat("lower_pad_",p_counter).c_str(),
			cat("lower_pad_",p_counter).c_str(),
			0,0,1,0.2);
  ++p_counter;
  p1_2->SetTopMargin(0.0025);
  p1_2->SetBottomMargin(0.3*bottom_marg_scale);
  p1_2->SetLeftMargin(0.1*can_size_ratio*left_marg_scale);
  p1_2->SetRightMargin(0.05);
  return {c,p1_1,p1_2};
}
DoubleCanvasPtr MakeDoubleCanvas(const char* title, int width, int height)
{
  static int c_counter = 0;
  TCanvas* c = new TCanvas(cat("DoubleCanvas_",c_counter++).c_str(),title,0,0,width,height);
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


void DrawDistribution(
		      CanvasPtr& canvas,
		      HistArray& histograms,
		      const HiggsModel& THDM,
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

      if (WRITE)
	{
	  histograms[0]->Write();
	  histograms[1]->Write();
	  if (NLO) histograms[3]->Write();
	}
      
      histograms[0]->DrawCopy("same hist ][");
      
      leg->AddEntry(histograms[0] ,"LO QCD","L");
      if (THDM.NBosons()>1)
	{
	  leg->AddEntry(histograms[1] ,"LO QCD + #Phi_{1} + #Phi_{2} ","L");
	}
      else
	{
	  leg->AddEntry(histograms[1] ,"LO QCD + #Phi ","L");
	}
      if (NLO)
	{
	  if (THDM.NBosons()>1)
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


void DrawTwoDistributions(
			  DoubleCanvasPtr& canvas,
			  HistArray& histogramsL,
			  HistArray& histogramsR,
			  const HiggsModel& THDM,
			  double* norm,
			  std::string const& titleLX,std::string const& titleRX,
			  std::string const& titleLY,std::string const& titleRY,
			  bool WRITE,
			  bool NLO)
{
  DrawDistribution(
		   canvas.c1,
		   histogramsL,
		   THDM,
		   norm,
		   titleLX,
		   titleLY,
		   WRITE,
		   NLO,
		   1);
  DrawDistribution(
		   canvas.c2,
		   histogramsR,
		   THDM,
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
  // hist->Smooth();
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
			 // top/antitop 4-momenta in the tt z.m.f.
			 FV const& k1,
			 FV const& k2,
			 // spin vectors in the tt z.m.f.
			 FV& s1,
			 FV& s2,
			 // spin vectors in t and tbar restframe, resp.
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


