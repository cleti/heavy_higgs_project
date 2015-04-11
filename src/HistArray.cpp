

#include "../inc/HistArray.h"






const char* cat(std::string str, int number)
{
  std::ostringstream strstr;
  strstr << str << number;
  return strstr.str().c_str();
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// histogram array counter to generate unique labels for each ROOT histogram
int HistArray::d_ID = 0;

HistArray::HistArray(int nbinsx,
		     double xlow,
		     double xup,
		     int mass_dim,
		     std::string const& label,
		     bool SUMW2):
  d_histograms({TH1D(cat("LO_QCD_",d_ID),"",nbinsx,xlow,xup),
	TH1D(cat("LO_INT_",d_ID),"",nbinsx,xlow,xup),
	TH1D(cat("LO_PHI_",d_ID),"",nbinsx,xlow,xup),
	TH1D(cat("NLO_V_",d_ID) ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_D_",d_ID) ,"",nbinsx,xlow,xup),
	TH1D(cat("NLO_R_",d_ID) ,"",nbinsx,xlow,xup)}),
  d_active(BOOST_BINARY(11 111 1)),
  d_active_t(0),
  d_label(label),
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

void HistArray::Normalize(const double& mScale)
{
  // relies on equidistant binning in all histograms !!!
  double binw = d_histograms[0].GetBinWidth(1);
  double fmass= std::pow(mScale,d_mass_dim);
  std::cout << "\n Normalizing: " << d_label;
  std::cout << "\n   binw: " << binw;
  std::cout << "\n   mass: " << d_mass_dim << " -> " << fmass << std::endl;  
  for (unsigned i=0;i<NHIST;i++)
    {
      d_histograms[i].Scale(fmass/(binw));
      std::cout << "  hist " << i <<  " done..\n";
    }
  std::cout << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////
double DIST_M_L = 340;
double DIST_M_U = 800;

double DIST_P_L = 60;
double DIST_P_U = 500;

double DIST_Y_L =-3.0;
double DIST_Y_U = 3.0;

double DIST_A_L = -1.0;
double DIST_A_U =  1.0;
////////////////////////////////

// prepare some histograms
HistArray MttDistributions(23,DIST_M_L,DIST_M_U,1,"M_{t#bar{t}} distribution");
HistArray PT1Distributions(22,DIST_P_L,DIST_P_U,1,"p_{T,t} distribution");
HistArray PT2Distributions(22,DIST_P_L,DIST_P_U,1,"p_{T,#bar{t}} distribution");
HistArray PT12Distributions(25,0.0,500.0,1,"p_{T,t+#bar{t}} distribution");
HistArray Y1Distributions(14,DIST_Y_L,DIST_Y_U,0,"#eta_{t} distribution");
HistArray Y2Distributions(14,DIST_Y_L,DIST_Y_U,0,"#eta_{#bar{t}} distribution");
HistArray DYDistributions(14,DIST_Y_L,DIST_Y_U,0,"#Delta |y| distribution");
// HistArray PTL1Distributions(40,DIST_P_L,DIST_P_U,1,"Lepton p_{T} distribution");
// HistArray PTL2Distributions(40,DIST_P_L,DIST_P_U,1,"Anti-Lepton p_{T} distribution");
// HistArray PHILLDistributions(20,DIST_A_L,DIST_A_U,0,"Lepton-Anti-lepton opening angle");
// HistArray SPartDistributions(40,DIST_M_L,DIST_M_U,1,"Partonic cross section");





/////////////////////////////////////////////////////////////////////////////////////////////
// define observables for differential distributions here ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// mass of 2-particle system [GeV]
double obs_M12(FV const& k1, FV const& k2)
{
  return sqrt(sp(k1,k1)+sp(k2,k2)+2.0*sp(k1,k2));
}
// transverse momentum (transverse plane is the x-y-plane) [GeV]
double obs_PT(FV const& k)
{
  return sqrt(k[1]*k[1]+k[2]*k[2]);
}
// pseudo-rapidity
double obs_Y(FV const& k)
{
  // the rapidity is not invariant under boosts in z-direction !
  // note that e.g. the reduced dipole phase spaces are boosted along z
  // with respect to the parton c.m.f.
  // double K  = LEN(k);
  // return 0.5*log(((K+k[3])/(K-k[3])));
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
  return (k1[1]*k2[1]+k1[2]*k2[2]+k1[3]*k2[3])/(K1*K2);
}
