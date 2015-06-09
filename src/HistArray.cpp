

#include "../inc/HistArray.h"



static std::string cat(std::string str, int number)
{
  std::ostringstream strstr;
  strstr << str << number;
  return strstr.str();
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// histogram array counter to generate unique labels for each ROOT histogram
int HistArray::d_ID = 0;
// layout settings for the generated TH1D's
double HistArray::LabelSizeY   = 0.028;
double HistArray::TitleOffsetY = 1.2;


HistArray::HistArray(int nbinsx,
		     double xlow,
		     double xup,
		     int mass_dim,
		     std::string const& name,
		     std::string const& lab_x,
		     std::string const& lab_y,
		     bool SUMW2):
  d_histograms({
        TH1D(cat("LO_QCD_", d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("NLO_QCD_",d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("LO_PHI_", d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_V_",d_ID).c_str() ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_D_",d_ID).c_str() ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_R_",d_ID).c_str() ,"",nbinsx,xlow,xup)}),
  d_active(0),
  d_active_t(0),
  d_name(name),
  d_label_x(),
  d_label_y(),
  d_mass_dim(mass_dim)

{
  for (int i=0;i<NHIST;++i)
    {
      d_histograms[i].GetYaxis()->SetLabelSize(HistArray::LabelSizeY);
      d_histograms[i].GetYaxis()->SetTitleOffset(HistArray::TitleOffsetY);
      d_histograms[i].GetXaxis()->SetTitle(lab_x.c_str());
      d_histograms[i].GetYaxis()->SetTitle(lab_y.c_str());
      d_histograms[i].SetLineWidth(1);
      d_histograms[i].SetLineColor(1);
      if (SUMW2)  d_histograms[i].Sumw2();
    }
  d_ID++;
}

void HistArray::Normalize(const double& mScale, int verb)
{
  double binc = 0.0; // temp
  double fmass= std::pow(mScale,d_mass_dim);
  if (verb>1)
    {
      std::cout << "\n Normalizing: " << d_name;
      std::cout << "\n   mass: " << d_mass_dim << " -> " << fmass << std::endl;
    }
  for (unsigned i=0;i<NHIST;i++)
    {
      // normalize only activated histograms
      if (IsActive(i))
	{
	  for (int j=1; j<=d_histograms[i].GetNbinsX(); ++j)
	    {

	      // ROOT counts bins from 1 to NBinsX, 0 and NBinsX+1 are underflow and overflow
	      binc = d_histograms[i].GetBinContent(j);
	      d_histograms[i].SetBinContent(j,binc*fmass/(d_histograms[i].GetBinWidth(j)));
	    }
	  if (verb>1) std::cout << "  hist " << i <<  " done..\n";
	}
      else
	{
	  if (verb>1) std::cout << "  hist " << i <<  " deactivated!\n";
	}
    }
  if (verb>1) std::cout << std::endl;
}



void HistArray::Print(std::ostream& ost)
{
  const int colWidth1 = 16;
  const int colWidth2 = 20;
  
  ost << std::endl;
  ost << "#  "  << d_name << std::endl;
  ost << "#  " << std::setw(colWidth1) << "Bin low edge";
  ost << std::setw(colWidth2) << "QCD [LO]";
  ost << std::setw(colWidth2) << "QCD+PHI [LO]";
  ost << std::setw(colWidth2) << "QCD [NLO]";
  ost << std::setw(colWidth2) << "QCD+PHI [NLO]" << std::endl << std::endl;

  // integrals
  double I_QCD_LO = 0.0;
  double I_QCD_PHI_LO = 0.0;
  double I_QCD_NLO = 0.0;
  double I_QCD_PHI_NLO = 0.0;

  // current bin content
  double QCD_LO = 0.0;
  double QCD_PHI_LO = 0.0;
  double QCD_NLO = 0.0;
  double QCD_PHI_NLO = 0.0;


  for (int i=0;i<=d_histograms[0].GetNbinsX()+1;++i)
    {
      // get current bin values
      QCD_LO      = d_histograms[H_LO_QCD].GetBinContent(i);
      QCD_PHI_LO = QCD_LO
	+ d_histograms[H_LO_PHI].GetBinContent(i);   
      QCD_NLO     = QCD_LO + d_histograms[H_NLO_QCD].GetBinContent(i);
      QCD_PHI_NLO = QCD_NLO
	+ d_histograms[H_LO_PHI].GetBinContent(i) 
	+ d_histograms[H_NLO_PHI_V].GetBinContent(i)
	+ d_histograms[H_NLO_PHI_ID].GetBinContent(i)
	+ d_histograms[H_NLO_PHI_R].GetBinContent(i);
      // update integral values (ignore overflow and underflow bin)
      if (i>0 && i<=d_histograms[0].GetNbinsX())
	{
	  double binw = d_histograms[H_LO_QCD].GetBinWidth(i);
	  I_QCD_LO      += QCD_LO*binw;
	  I_QCD_PHI_LO  += QCD_PHI_LO*binw;
	  I_QCD_NLO     += QCD_NLO*binw;
	  I_QCD_PHI_NLO += QCD_PHI_NLO*binw;
	}
      // print bin values
      if (i>0)
	{
	  ost << std::setw(colWidth2) << std::setprecision(5)  << d_histograms[0].GetBinLowEdge(i);
	}
      else
	{
	  ost << std::setw(colWidth2) << std::setprecision(5)  << "";
	}
      ost << std::setw(colWidth2) << std::setprecision(10) << QCD_LO;
      ost << std::setw(colWidth2) << std::setprecision(10) << QCD_PHI_LO;
      ost << std::setw(colWidth2) << std::setprecision(10) << QCD_NLO;
      ost << std::setw(colWidth2) << std::setprecision(10) << QCD_PHI_NLO << std::endl;
    }
  
  ost << std::endl;
  ost << "#  " << std::setw(colWidth1) << "Integral: ";
  ost << std::setw(colWidth2) << I_QCD_LO;
  ost << std::setw(colWidth2) << I_QCD_PHI_LO;
  ost << std::setw(colWidth2) << I_QCD_NLO;
  ost << std::setw(colWidth2) << I_QCD_PHI_NLO << std::endl << std::endl;
}


void HistArray::Status(std::ostream& ost)
{
  ost << std::endl;
  ost << d_name << std::endl;
  for (int i=0;i<NHIST;++i)
    {
      ost << "#" << std::setw(2) << i << ": " << (IsActive(i)?'1':'0') << std::endl;
    }
  ost << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////
double DIST_M_L  = 340;
double DIST_M_U  = 1000;
double DIST_M_U2 = 900;

double DIST_P_L = 60;
double DIST_P_U = 500;

double DIST_Y_L = -3.0;
double DIST_Y_U =  3.0;

double DIST_A_L = -1.0;
double DIST_A_U =  1.0;
////////////////////////////////

// some predefined histograms
// spin-INdependent observables
HistArray Mtt_Histograms(33,DIST_M_L,DIST_M_U,1,
			 "Top/Antitop invariant mass distribution",
			 "M_{t#bar{t}} [GeV]",
			 "#frac{d#sigma}{dM_{t#bar{t}}} [pb/GeV]");

HistArray PT1_Histograms(22,DIST_P_L,DIST_P_U,1,
			 "Top transverse momentum distribution",
			 "p_{T,t} [GeV]",
			 "#frac{d#sigma}{dp_{T,t}} [pb/GeV]");
HistArray PT2_Histograms(22,DIST_P_L,DIST_P_U,1,
			 "Antitop transverse momentum distribution",
			 "p_{T,#bar{t}} [GeV]",
			 "#frac{d#sigma}{dp_{T,#bar{t}}} [pb/GeV]");
HistArray PT12_Histograms(10,0.0,500.0,1,
			  "Top+Antitop transverse momentum distribution",
			  "|p_{T,t}+p_{T,#bar{t}}| [GeV]",
			  "#frac{d#sigma}{d|p_{T,t}+p_{T,#bar{t}}|} [pb/GeV]");
HistArray Y1_Histograms(12,DIST_Y_L,DIST_Y_U,0,
			"Top rapidity distribution",
			"y_{t}",
			"#frac{d#sigma}{dy_{t}} [pb]");
HistArray Y2_Histograms(12,DIST_Y_L,DIST_Y_U,0,
			"Antitop rapidity distribution",
			"y_{#bar{t}}",
			"#frac{d#sigma}{dy_{#bar{t}}} [pb]");
HistArray DY_Histograms(12,DIST_Y_L,DIST_Y_U,0,
			"Distribution of top/antitop rapidity difference",
			"#Delta |y|",
			"#frac{d#sigma}{d #Delta |y|} [pb]");

const int NbinsS = 56;
// spin-dependent observables
HistArray PHIT12_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			    "Lepton-Anti-lepton transversal opening angle [lab frame]",
			    "M_{t#bar{t}} [GeV]",
			    "D_{T,open}^{lab} [pb/GeV]",
			    true);
HistArray Dopen_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			   "Lepton-Anti-lepton opening angle [ttbar z.m.f.]",
			   "M_{t#bar{t}} [GeV]",
			   "D_{open} [pb/GeV]",
			   true);
HistArray OCP1_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			  "CP-odd triple correlation",
			  "M_{t#bar{t}} [GeV]",
			  "O_{CP1} [pb/GeV]",
			  true);
HistArray B1_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			"Longitudinal polarization",
			"M_{t#bar{t}} [GeV]",
			"B_{1} [pb/GeV]",
			true);
HistArray Chel_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			  "Helicity angle distribution",
			  "M_{t#bar{t}} [GeV]",
			  "C_{hel} [pb/GeV]",
			  true);



