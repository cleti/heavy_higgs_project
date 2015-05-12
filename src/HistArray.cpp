

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

HistArray::HistArray(int nbinsx,
		     double xlow,
		     double xup,
		     int mass_dim,
		     std::string const& name,
		     bool SUMW2):
  d_histograms({TH1D(cat("LO_QCD_",d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("NLO_QCD_",d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("LO_PHI_",d_ID).c_str(),"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_V_",d_ID).c_str() ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_D_",d_ID).c_str() ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_PHI_R_",d_ID).c_str() ,"",nbinsx,xlow,xup)}),
  d_active(BOOST_BINARY(11 111 1)),
  d_active_t(0),
  d_name(name),
  d_label_x(),
  d_label_y(),
  d_mass_dim(mass_dim),
  d_obs(nullptr)
{
  if (SUMW2)
    {
      for (int i=0;i<NHIST;++i)
	{
	  d_histograms[i].Sumw2();
	}
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
  ost << std::endl;
  ost << "#  "  << d_name << std::endl;
  ost << "#  " << std::setw(17) << "Bin low edge" << std::setw(22) << "QCD [LO]" << std::setw(22) << "QCD+PHI [LO]" << std::setw(22) << "QCD+PHI [NLO]" << std::endl << std::endl;

  // integrals
  double I_QCD = 0.0;
  double I_QCD_PHI_LO = 0.0;
  double I_QCD_PHI_NLO = 0.0;

  // current bin content
  double QCD = 0.0;
  double QCD_PHI_LO = 0.0;
  double QCD_PHI_NLO = 0.0;

  // underflow bin
  QCD         = d_histograms[0].GetBinContent(0);
  QCD_PHI_LO  = QCD + d_histograms[1].GetBinContent(0);
  QCD_PHI_NLO = QCD_PHI_LO
    + d_histograms[3].GetBinContent(0)
    + d_histograms[4].GetBinContent(0)
    + d_histograms[5].GetBinContent(0);
  ost << std::setw(20) << std::setprecision(5)  << "" << "  ";
  ost << std::setw(20) << std::setprecision(10) << QCD    << "  ";
  ost << std::setw(20) << std::setprecision(10) << QCD_PHI_LO << "  ";
  ost << std::setw(20) << std::setprecision(10) << QCD_PHI_NLO << std::endl;
  for (int i=1;i<=d_histograms[0].GetNbinsX()+1;++i)
    {
      // get current bin values
      QCD         = d_histograms[0].GetBinContent(i);
      QCD_PHI_LO  = QCD + d_histograms[1].GetBinContent(i);
      QCD_PHI_NLO = QCD_PHI_LO
	+ d_histograms[3].GetBinContent(i)
	+ d_histograms[4].GetBinContent(i)
	+ d_histograms[5].GetBinContent(i);
      // update integral values (ignore overflow bin)
      if (i<=d_histograms[0].GetNbinsX())
	{
	  I_QCD         += QCD*d_histograms[0].GetBinWidth(i);
	  I_QCD_PHI_LO  += QCD_PHI_LO*d_histograms[1].GetBinWidth(i);
	  I_QCD_PHI_NLO += QCD_PHI_NLO*d_histograms[3].GetBinWidth(i);
	}
      // print bin values
      ost << std::setw(20) << std::setprecision(5)  << d_histograms[0].GetBinLowEdge(i) << "  ";
      ost << std::setw(20) << std::setprecision(10) << QCD    << "  ";
      ost << std::setw(20) << std::setprecision(10) << QCD_PHI_LO << "  ";
      ost << std::setw(20) << std::setprecision(10) << QCD_PHI_NLO    << std::endl;
    }
  ost << std::endl;
  ost << "#  " << std::setw(17) << "Integral: " << std::setw(22) << I_QCD << std::setw(22) << I_QCD_PHI_LO << std::setw(22) << I_QCD_PHI_NLO << std::endl << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////
double DIST_M_L = 340;
double DIST_M_U = 1000;

double DIST_P_L = 60;
double DIST_P_U = 500;

double DIST_Y_L = -3.0;
double DIST_Y_U =  3.0;

double DIST_A_L = -1.0;
double DIST_A_U =  1.0;
////////////////////////////////

// prepare some histograms
// spin independent observables
HistArray MttDistributions(33,DIST_M_L,DIST_M_U,1,"M_{t#bar{t}} distribution");
HistArray PT1Distributions(22,DIST_P_L,DIST_P_U,1,"p_{T,t} distribution");
HistArray PT2Distributions(22,DIST_P_L,DIST_P_U,1,"p_{T,#bar{t}} distribution");
HistArray PT12Distributions(10,0.0,500.0,1,"p_{T,t+#bar{t}} distribution");
HistArray Y1Distributions(14,DIST_Y_L,DIST_Y_U,0,"#eta_{t} distribution");
HistArray Y2Distributions(14,DIST_Y_L,DIST_Y_U,0,"#eta_{#bar{t}} distribution");
HistArray DYDistributions(14,DIST_Y_L,DIST_Y_U,0,"#Delta |y| distribution");

// angular distributions (spin dependent)
// HistArray PHIT12Distributions(12,-1.0,1.0,0,"Lepton-Anti-lepton transversal opening angle [lab frame]",true);
// HistArray PHI12Distributions(12,-1.0,1.0,0,"Lepton-Anti-lepton opening angle",true);
// HistArray OCPDistributions(12,-1.0,1.0,1,"CP-odd triple correlation",true);
// HistArray LP1Distributions(12,-1.0,1.0,1,"Top longitudinal polarization",true);
// HistArray LP2Distributions(12,-1.0,1.0,1,"Antitop longitudinal polarization",true);
// HistArray OHELDistributions(12,-1.0,1.0,1,"Helicity angle correlation",true);
HistArray PHIT12Distributions(33,DIST_M_L,DIST_M_U,1,"Lepton-Anti-lepton transversal opening angle [lab frame]",true);
HistArray DopenDistributions(33,DIST_M_L,DIST_M_U,1,"Lepton-Anti-lepton opening angle",true);
HistArray OCPDistributions(33,DIST_M_L,DIST_M_U,1,"CP-odd triple correlation",true);
HistArray ChelDistributions(33,DIST_M_L,DIST_M_U,1,"Top/Antitop longitudinal polarization",true);




/////////////////////////////////////////////////////////////////////////////////////////////
// define observables for differential distributions here ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// mass of 2-particle system [GeV]
double obs_M12(FV const& k1, FV const& k2)
{
  return sqrt(sp(k1,k1)+sp(k2,k2)+2.0*sp(k1,k2));
}
// transverse momentum (transverse plane is the x-y-plane) 
double obs_PT(FV const& k)
{
  return sqrt(k[1]*k[1]+k[2]*k[2]);
}

// transverse momentum of two particle system |k_{T,1}+k_{T,2}| 
double obs_PT12(FV const& k1, FV const& k2)
{
  return sqrt(pow(k1[1]+k2[1],2)+pow(k1[2]+k2[2],2)); 
}
// pseudo-rapidity
double obs_Y(FV const& k)
{
  return 0.5*log(((k[0]+k[3])/(k[0]-k[3])));
}
// difference of abs. values of transverse rapidities
double obs_DY(FV const& k1, FV const& k2)
{
  return fabs(obs_Y(k1))-fabs(obs_Y(k2));
}
// 3-dim opening angle of the two vectors
double obs_PHI(FV const& k1, FV const& k2)
{
  double K1 = LEN(k1);
  double K2 = LEN(k2);
  // PRINT(K1);
  // PRINT_4VEC(1.0/K1*k1);
  // PRINT(K2);
  // PRINT_4VEC(1.0/K2*k2);
  return (k1[1]*k2[1]+k1[2]*k2[2]+k1[3]*k2[3])/(K1*K2);
}
// 2-dim transversal opening angle of the two vectors
double obs_PHIT(FV const& k1, FV const& k2)
{
  double K1 = sqrt(k1[1]*k1[1]+k1[2]*k1[2]);
  double K2 = sqrt(k2[1]*k2[1]+k2[2]*k2[2]);
  return (k1[1]*k2[1]+k1[2]*k2[2])/(K1*K2);
}
// spatial triple product of three 4-vectors k1,k2,k3
double obs_TriProd(FV const& k1, FV const& k2, FV const& k3)
{
  // = (k1 x k2) dot k3 / |k1 x k2| / |k3|
  double t3 = k1[0] * k2[1] - k1[1] * k2[0];
  double t7 = -k1[0] * k2[2] + k1[2] * k2[0];
  double t11 = k1[1] * k2[2] - k1[2] * k2[1];
  double t14 = t3 * t3;
  double t15 = t7 * t7;
  double t16 = t11 * t11;
  double t18 = std::pow(k3[0], 2);
  double t19 = std::pow(k3[1], 2);
  double t20 = std::pow(k3[2], 2);
  double t23 = std::sqrt((t14 + t15 + t16) * (t18 + t19 + t20));
  return (t11 * k3[0] + t3 * k3[2] + t7 * k3[1]) / t23;
}

