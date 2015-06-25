

#include "../inc/HistArray.h"




const H_IndexMap imap_default = {
  {0,H_LO_QCD},
  {1,H_NLO_QCD},
  {2,H_LO_PHI},
  {3,H_NLO_PHI_V},
  {4,H_NLO_PHI_ID},
  {5,H_NLO_PHI_R},
};


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
      //d_histograms[i].GetYaxis()->SetMaxDigits(3);
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
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_LO_QCD].GetBinContent(i);
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_NLO_QCD].GetBinContent(i);
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_LO_PHI].GetBinContent(i);
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_NLO_PHI_V].GetBinContent(i);  
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_NLO_PHI_ID].GetBinContent(i);  
      ost << std::setw(colWidth2) << std::setprecision(10) << d_histograms[H_NLO_PHI_R].GetBinContent(i) << std::endl;    
      // ost << std::setw(colWidth2) << std::setprecision(10) << QCD_LO;
      // ost << std::setw(colWidth2) << std::setprecision(10) << QCD_PHI_LO;
      // ost << std::setw(colWidth2) << std::setprecision(10) << QCD_NLO;
      // ost << std::setw(colWidth2) << std::setprecision(10) << QCD_PHI_NLO << std::endl;
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



int HistArray::ReadFile(H_Index I, const boost_path& path, int tabNum, int colNum)
{
  const int colWidth1 = 16;
  const int colWidth2 = 20;
  PRINT(I);
  
  // if (!file) std::make_shared<std::ifstream>(path.string())
  std::ifstream file(path.string());
  while (!file.is_open())
    {
      FileBrowser fb(path.parent_path().string());
      std::cout << std::endl << " Could not open file " << path << std::endl;
      std::cout << " Pick new input file:" << std::endl;
      file.open(fb.browse());
    }

  // alias for the current histgram in 'DistVec* dist_vec'
  TH1D& hist  = d_histograms[I];

  std::string sbuff;
  int rowCount = 0;
  int tabCount = 0;
  int rowMax   = hist.GetNbinsX();
  double bin;
  double val;
  while(file.good())
    {
      if (rowCount>rowMax+1)
	{
	  std::cout << "\n\n Error! ";
	  std::cout << " Table dimension of input data does not match the histogram!\n\n ";
	  return 0;
	}
      std::getline(file,sbuff);
      // remove leading whitespace
      while (sbuff.front() == ' ' || sbuff.front() == '\t')
	{
	  sbuff.erase(std::begin(sbuff));
	}
      // inspect line
      if (sbuff.size()>0) 
	{
	  if (sbuff.front() == '#')
	    { // comment in combination with rowCount>0 signals next distribution:
	      // cout nbins and reset, print comment, increase tabCount
	      if (rowCount>0)
		{
		  std::cout << std::endl;
		  std::cout << " [ END OF TABLE " << tabCount;
		  std::cout << ", " << rowCount << " bins ";
		  if (tabCount != tabNum)
		    {
		      std::cout << " ignored ";
		    }
		  std::cout << "]\n" << std::endl;
		  rowCount = 0;
		  ++tabCount;
		}
	      if (tabCount == (tabNum-1)) std::cout << sbuff << std::endl;
	    }
	  else
	    { // not a comment: extract numbers, start counting bins if this is the right tableNum
	      ++rowCount;
	      if (tabCount == tabNum)
		{
		  // if there is no space or something after the last number
		  // iss.good() will be false and the last number wont be copied
		  sbuff+=' ';
		  std::stringstream iss(sbuff);
		  ///////////////////////////////
		  // first column is lower bin boundary
		  iss >> bin;
		  std::cout << std::setw(5) << rowCount << " ";
		  std::cout << std::setw(colWidth1) << bin << " ";	      
		  int colCount = 0;
		  // extract the subsequent columns
		  while (1)
		    {
		      iss >> val;
		      // copy value if this is the right colNum
		      // normalize to bin width
		      
		      if (iss.good())
			{
			  if (colCount == colNum)
			    {
			      double binw = hist.GetBinWidth(rowCount);
			      hist.SetBinContent(rowCount,val/binw*4.0/81.0);
			    }
			  std::cout << std::setw(colWidth2) << val << " ";
			  ++colCount;
			}
		      else
			{
			  break;
			}
		    }
		  std::cout << std::endl;
		  ///////////////////////////////

		  if ( std::fabs( bin - (hist.GetBinLowEdge(rowCount)) ) > 1e-5 )
		    {
		      std::cout << "\n\n Error at bin (low-edge hist): " << hist.GetBinLowEdge(rowCount);
		      std::cout << "\n\n Error at bin (low-edge input):" << bin;
		      std::cout << "\n Input data does not match histogram binning!\n\n";
		      // discard data for this histogram
		      hist.Reset();
		      break;
		    }
		}
	    }
	}
    }
  if (file.eof())
    {
      if (rowCount>0)
	{
	  std::cout << std::endl;
	  std::cout << " [ END OF TABLE " << tabCount;
	  std::cout << ", " << rowCount << " bins ";
	  if (tabCount != tabNum)
	    {
	      std::cout << " ignored ";
	    }
	  std::cout << "]\n" << std::endl;
	  rowCount = 0;
	  ++tabCount;
	}     
      std::cout << "\n End-of-File reached on input operation.\n " << std::endl;
      return 1;
    }
  if (file.fail())
    {
      std::cout << std::endl << " Error! ";
      std::cout << "\n Logical error on i/o operation.\n " << std::endl;
      return 0;
    }
  if (file.bad())
    {
      std::cout << std::endl << " Error! ";
      std::cout << "\n Read/writing error on i/o operation.\n " << std::endl;
      return 0;
    }
  return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
int HistArray::ReadTable(
			 std::ifstream& file,
			 const double& norm,
			 const H_IndexMap& imap,
			 bool discardCurrentTable)
{
  if (!file.is_open())
    {
      return 0;
    }
  // const int colWidth1 = 16;
  // const int colWidth2 = 20;
  // alias for the current histgram in 'DistVec* dist_vec'
  TH1D&  hist0 = d_histograms[H_LO_QCD];
  std::string buff;
  //std::stringstream sbuff;
  int rowCount = 0;
  int rowMax   = hist0.GetNbinsX();
  double binLowEdge;
  double binVal;

  // MAIN LOOP ////////////////////////////////////////////////////////////////////////////////
  while(file.good())
    {
      // FETCH INPUT LINE /////////////////////////////////////////////////////////////////////
      std::getline(file,buff);
      // if there is no space or something after the last number
      // iss.good() will be false already after extraction of the last number
      buff+=' ';
      // remove leading whitespace
      while (buff.front() == ' ' || buff.front() == '\t')
	{
	  buff.erase(std::begin(buff));
	}            
      // INSPECT LINE ////////////////////////////////////////////////////////////////////////
      if (buff.size()>0) 
	{
	  if (buff.front() == '#')
	    { // comment-line: comment in combination with rowCount>0 signals end-of-table
	      if (rowCount>0)
		{
		  std::cout << std::endl;
		  std::cout << " [ END OF TABLE ";
		  std::cout << ", " << rowCount << " rows ]\n\n";
		  return 1;
		}
	      if (discardCurrentTable) discardCurrentTable = false;
	    }
	  else
	    { // data-line: extract numbers, count table row and compare with target histogram size
	      // discardCurrentTable: do nothing until a comment-line was processed
	      if (discardCurrentTable) continue;
	      // use stringstream to extract numbers
	      std::stringstream sbuff(buff);
	      // first column is the lower bin boundary
	      sbuff >> binLowEdge;    
	      ////////////////////////////////////////////////////////////////////////////////////
	      if ( ++rowCount > rowMax+1)
		{
		  std::cout << "\n Error reading row " << rowCount;
		  std::cout << "\n  Buffer content (str): " << buff;
		  std::cout << "\n  Too many rows in input table!\n\n ";
		  return 0;
		}	      
	      // compare bin boundary given in current row with target histogram binning
	      if ( std::fabs( binLowEdge - (hist0.GetBinLowEdge(rowCount)) ) > 1e-5 )
		{
		  std::cout << "\n Error reading row " << rowCount;
		  std::cout << "\n  Buffer content (str): " << buff;
		  std::cout << "\n  low-edge hist: " << hist0.GetBinLowEdge(rowCount);
		  std::cout << "\n  low-edge input:" << binLowEdge;
		  std::cout << "\n Input data does not match histogram binning!\n\n";
		  // discard data from all histograms
		  //this->Reset();
		  return 0;
		}	      
	      ////////////////////////////////////////////////////////////////////////////////////
	      // std::cout << std::setw(5) << rowCount << " ";
	      // std::cout << std::setw(colWidth1) << binLowEdge << " ";	  
	      int colCount = 0;
	      // extract the subsequent columns
	      while (1)
		{
		  sbuff >> binVal;
		  if (sbuff.good())
		    {
		      // std::cout << std::setw(colWidth2) << binVal << " ";
		      // get histogram index corresponding to colCount
		      H_IndexMap::const_iterator it = imap.find(colCount++);
		      // fill value into histogram if the index was found in the map
		      if (it != imap.end())
			{
			  // std::cout << std::setw(colWidth2/2) << it->second << " ";
			  // alias for the current histgram in 'DistVec* dist_vec'
			  TH1D&  hist = d_histograms[it->second];
			  // zongguos numbers are not normalized to bin width and
			  // on tt bar level -> factor 4/81
			  hist.SetBinContent(rowCount,binVal/hist.GetBinWidth(rowCount)*norm);
			}
		    }
		  else
		    {
		      break;
		    }
		}
	      // std::cout << std::endl;  
	    }
	}
      // !INSPECT LINE ////////////////////////////////////////////////////////////////////////
    }
  // !MAIN LOOP ///////////////////////////////////////////////////////////////////////////////
  if (file.eof())
    {
      if (rowCount>0)
	{
	  std::cout << std::endl;
	  std::cout << " [ END OF TABLE ";
	  std::cout << ", " << rowCount << " rows ]\n\n";

	}     
      std::cout << "\n End-of-File reached on input operation.\n " << std::endl;
      return 0;
    }
  if (file.fail())
    {
      std::cout << "\n Logical error on i/o operation.";
      return 0;
    }
  if (file.bad())
    {
      std::cout << "\n Read/writing error on i/o operation.\n\n";
      return 0;
    }
  return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////
double DIST_M_L  = 340;
double DIST_M_U  = 1000;
double DIST_M_U2 = 1200;

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

const int NbinsS = 86;
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
			"Top longitudinal polarization",
			"M_{t#bar{t}} [GeV]",
			"B_{1} [pb/GeV]",
			true);
HistArray B2_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			"Antitop longitudinal polarization",
			"M_{t#bar{t}} [GeV]",
			"B_{2} [pb/GeV]",
			true);
HistArray Chel_Histograms(NbinsS,DIST_M_L,DIST_M_U2,1,
			  "Helicity angle distribution",
			  "M_{t#bar{t}} [GeV]",
			  "C_{hel} [pb/GeV]",
			  true);






int ReadData(
	     std::string path,
	     DistVec::iterator it,
	     DistVec::const_iterator it_end,
	     const double& norm,
	     const H_IndexMap& imap)
{
  // open file with input data
  std::ifstream file(path);
  while (!file.is_open())
    {
      FileBrowser fb(path);
      std::cout << std::endl << " Could not open file " << path << std::endl;
      std::cout << " Pick new input file:" << std::endl;
      file.open(path=fb.browse());
    }
  std::cout << "\n\n Reading data from file '" << path << "' ...\n\n";
  
  int ret = 1;
  int tabCount = 0;
  // the tables in the data file should be in the same order as the distributions in distVec!
  for ( ; it != it_end; ++it)
    {
      // skip first distribution -> Mtt
      // if (dist==distVec.begin()) continue;
      std::cout << " Processing table " << (tabCount++) << " -> '" << (*it)->GetName() << "' ...\n";
      ret *= (*it)->GetHistograms()->ReadTable(file,norm,imap,!ret);
    }
  std::cout << std::endl;
  return ret;
}




