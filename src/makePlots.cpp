



#include "../inc/Makros.h"

#include "boost/filesystem.hpp" 
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
using namespace boost::filesystem;    
using namespace std;

// ROOT header files
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>



// directory layout for each result data set:
const string scale_m_path = "mu1.0";
const string scale_h_path = "mu2.0";
const string scale_l_path = "mu0.5";
// each of these should contain a ROOT file with the histograms to be plotted

// some constants to compute the event rates in ll and lj channels
const double eps_ll = 0.22;
const double Br_ll  = 4.0/81.0;
const double eps_lj = 0.12;
const double Br_lj  = 24.0/81.0;
double Lumi   = 100*1000.0; // integrated lumi of 8TeV run in [pb^-1]


TH1D* makeErrorHist(TH1D* hist_l, TH1D* hist_h)
{
  if (!hist_l || !hist_h)
    {
      cout << endl << " In " << __FUNCTION__ << ": nullpointer received." << endl;
      exit(0);
    }
  if (hist_l->GetNbinsX() != hist_h->GetNbinsX())
    {
      cout << endl << " In " << __FUNCTION__ << ": binning of input histograms does not match." << endl;
      exit(0);
    }

  TH1D* ret = new TH1D(*hist_l);
  ret->SetNameTitle("scale_band","scale_band");
  for (int i=1;i<=hist_l->GetNbinsX();++i)
    {
      double cnt = (hist_l->GetBinContent(i)+hist_h->GetBinContent(i))/2.0;
      double err = 2.0*fabs(hist_h->GetBinContent(i)-cnt);
      ret->SetBinContent(i,cnt);
      ret->SetBinError(i,err);
    }

  ret->SetFillStyle(3002);
  ret->SetFillColor(1);
  ret->SetLineWidth(0);
  return ret;  
}








int main(int argc, char** argv)
{
  string results_path_name = "results/";
  if (argc > 1)
    {
      results_path_name = argv[1];
    }

  path results(results_path_name);

  // request valid directory containting results
  while (!exists(results))
    {
      cout << endl << " Directory '" << results_path_name;
      cout << "' does not exist. Enter new path: " << endl; 
      cout << " >> "; cin >> results_path_name;
      results = results_path_name;
    }

  // choose data set to be plotted
  int user_sel = -1;
  vector<directory_entry> dir_entries;
  directory_iterator end_itr; // default construction yields past-the-end
  cout << endl << " Directories in " << results << ":" << endl;
  while (user_sel < 0)
    {
      int cnt = 0;
      for ( directory_iterator itr(results);
	    itr != end_itr;
	    ++itr )
	{
	  if ( is_directory(itr->status()) )
	    {
	      cout << setw(5) << cnt << setw(20) << itr->path().filename() << endl;
	      dir_entries.push_back(*itr);
	      ++cnt;
	    }
	}
      cout << endl << " Choose directory: " << endl; 
      cout << " >> "; cin >> user_sel;
      if (user_sel>=cnt) user_sel = -1;
    }

  // get the ROOT files
  cout << endl <<  " User selection: " << dir_entries[user_sel].path() << endl;


  
  // open .root file with histograms mu=1
  string filename_scale_m;
  path path_tmp = dir_entries[user_sel].path()/scale_m_path;
  if (!exists(path_tmp))
    {
      cout << endl << " Path: " << path_tmp << " does not exist." << endl;
      exit(0);
    }

  for ( directory_iterator itr(path_tmp);
	itr != end_itr;
	++itr )
    {
      // get first .root file
      if(is_regular_file(*itr) && itr->path().extension() == ".root")
	{
	  filename_scale_m = itr->path().string();
	  break;
	}
    }
  if (filename_scale_m.empty())
    {
      cout << endl << " No .root file found in : " << path_tmp << "." << endl;
      exit(0);
    }
  
  // open .root file with histograms mu=2
  string filename_scale_h;
  path_tmp = dir_entries[user_sel].path()/scale_h_path;
  if (!exists(path_tmp))
    {
      cout << endl << " Path: " << path_tmp << " does not exist." << endl;
      exit(0);
    }

  for ( directory_iterator itr(path_tmp);
	itr != end_itr;
	++itr )
    {
      // get first .root file
      if(is_regular_file(*itr) && itr->path().extension() == ".root")
	{
	  filename_scale_h = itr->path().string();
	  break;
	}
    }
  if (filename_scale_h.empty())
    {
      cout << endl << " No .root file found in : " << path_tmp << "." << endl;
      exit(0);
    }
    
  // open .root file with histograms mu=0.5
  string filename_scale_l;
  path_tmp = dir_entries[user_sel].path()/scale_l_path;
  if (!exists(path_tmp))
    {
      cout << endl << " Path: " << path_tmp << " does not exist." << endl;
      exit(0);
    }

  for ( directory_iterator itr(path_tmp);
	itr != end_itr;
	++itr )
    {
      // get first .root file
      if(is_regular_file(*itr) && itr->path().extension() == ".root")
	{
	  filename_scale_l = itr->path().string();
	  break;
	}
    }
  if (filename_scale_l.empty())
    {
      cout << endl << " No .root file found in : " << path_tmp << "." << endl;
      exit(0);
    }


  cout << endl << " mu = 1.0 data:" << endl << filename_scale_m;
  cout << endl << " mu = 2.0 data:" << endl << filename_scale_h;
  cout << endl << " mu = 0.5 data:" << endl << filename_scale_l;
  cout << endl;

  TFile* file_scale_m = TFile::Open(filename_scale_m.c_str(),"READ");
  TFile* file_scale_h = TFile::Open(filename_scale_h.c_str(),"READ");
  TFile* file_scale_l = TFile::Open(filename_scale_l.c_str(),"READ");

  file_scale_m->ls();
  file_scale_h->ls();
  file_scale_l->ls();

  TH1D* h_lo_qcd_scale_m  = nullptr;
  TH1D* h_lo_qcd_scale_h  = nullptr;
  TH1D* h_lo_qcd_scale_l  = nullptr;
      
  TH1D* h_lo_phi_scale_m  = nullptr;
  TH1D* h_lo_phi_scale_h  = nullptr;
  TH1D* h_lo_phi_scale_l  = nullptr;
      
  TH1D* h_nlo_qcd_scale_m = nullptr;
  TH1D* h_nlo_qcd_scale_h = nullptr;
  TH1D* h_nlo_qcd_scale_l = nullptr;
      
  TH1D* h_nlo_phi_scale_m = nullptr;
  TH1D* h_nlo_phi_scale_h = nullptr;
  TH1D* h_nlo_phi_scale_l = nullptr;
  
  // flags indicating if NLO histograms are available
  bool NLO_QCD = true;
  bool NLO_PHI = true;
      
  int obj_num = -1;
  string histName_qcd_base = "LO_QCD_";
  string histName_phi_base = "LO_PHI_";
  while (true)
    {
      cout << endl << " Choose distribution: " << endl; 
      cout << " >> "; cin >> obj_num;
   
      stringstream histName_lo_qcd; histName_lo_qcd << histName_qcd_base << obj_num;
      stringstream histName_lo_phi; histName_lo_phi << histName_phi_base << obj_num;
  
      // retrieve LO histograms from the .root files

      file_scale_m->GetObject(histName_lo_qcd.str().c_str(),h_lo_qcd_scale_m);
      file_scale_h->GetObject(histName_lo_qcd.str().c_str(),h_lo_qcd_scale_h);
      file_scale_l->GetObject(histName_lo_qcd.str().c_str(),h_lo_qcd_scale_l);
  
      file_scale_m->GetObject(histName_lo_phi.str().c_str(),h_lo_phi_scale_m);
      file_scale_h->GetObject(histName_lo_phi.str().c_str(),h_lo_phi_scale_h);
      file_scale_l->GetObject(histName_lo_phi.str().c_str(),h_lo_phi_scale_l);

      // LO histograms are required
      if (!h_lo_qcd_scale_m ||
	  !h_lo_qcd_scale_h ||
	  !h_lo_qcd_scale_l)
	{
	  cout << endl << " Error reading " << histName_qcd_base << " histograms from .root files. Check layout. " << endl;
	  histName_qcd_base = "NLO_QCD_";
	  continue;
	}
      if (!h_lo_phi_scale_m ||
	  !h_lo_phi_scale_h ||
	  !h_lo_phi_scale_l)
	{
	  cout << endl << " Error reading " << histName_phi_base << " histograms from .root files. Check layout. " << endl;
	  histName_phi_base = "NLO_PHI_";
	  continue;
	}


  
      stringstream histName_nlo_qcd; histName_nlo_qcd << "NLO_QCD_"   << obj_num;
      stringstream histName_nlo_phi; histName_nlo_phi << "NLO_PHI_V_" << obj_num;

      // retrieve NLO histograms from the .root files  
      file_scale_m->GetObject(histName_nlo_qcd.str().c_str(),h_nlo_qcd_scale_m);  
      file_scale_h->GetObject(histName_nlo_qcd.str().c_str(),h_nlo_qcd_scale_h);
      file_scale_l->GetObject(histName_nlo_qcd.str().c_str(),h_nlo_qcd_scale_l);
  
      file_scale_m->GetObject(histName_nlo_phi.str().c_str(),h_nlo_phi_scale_m);  
      file_scale_h->GetObject(histName_nlo_phi.str().c_str(),h_nlo_phi_scale_h);
      file_scale_l->GetObject(histName_nlo_phi.str().c_str(),h_nlo_phi_scale_l);
  
      // NLO histograms are optional
      if (!h_nlo_qcd_scale_m ||
	  !h_nlo_qcd_scale_h ||
	  !h_nlo_qcd_scale_l)
	{
	  cout << endl << " No NLO QCD histograms available. " << endl;
	  NLO_QCD = false;
	}
      
      if (!h_nlo_phi_scale_m ||
	  !h_nlo_phi_scale_h ||
	  !h_nlo_phi_scale_l)
	{
	  cout << endl << " No NLO 2HDM histograms available. " << endl;
	  NLO_PHI = false;
	}
      break;
    }

  NLO_PHI = false;
  
  double xlo = 0.0;
  double xhi = 0.0;
  double xmin= (NLO_QCD)?h_nlo_qcd_scale_m->GetXaxis()->GetXmin():h_lo_qcd_scale_m->GetXaxis()->GetXmin();
  double xmax =(NLO_QCD)?h_nlo_qcd_scale_m->GetXaxis()->GetXmax():h_lo_qcd_scale_m->GetXaxis()->GetXmax();
  bool subIntervals = true;
  while (true)
    {
      cout << endl << " Enter x-range [x_min=";
      cout << xmin;
      cout << ",x_max=" << xmax << "], q/Q to skip: " << endl; 
      cout << " >> ";
      char c = 'a';
      cin.get(c);
      if (c=='q' || c=='Q')
	{
	  break;
	}
      else if (c=='I')
	{
	  subIntervals = false;
	}
      else if (c=='i')
	{
	  subIntervals = true;
	}
      else
	{
	  cout << endl << " Undefined input command. " << endl;
	  continue;
	}
      //cin.putback(c);
      cin >> xlo;
      cin >> xhi;
      
      if (xlo<xhi && xlo>=xmin && xhi<=xmax)
  	{
  	  int bin_lo = (NLO_QCD)?h_nlo_qcd_scale_m->GetXaxis()->FindBin(xlo):h_lo_qcd_scale_m->GetXaxis()->FindBin(xlo);
  	  int bin_hi = (NLO_QCD)?h_nlo_qcd_scale_m->GetXaxis()->FindBin(xhi):h_lo_qcd_scale_m->GetXaxis()->FindBin(xhi);
  	  bool ratioGtOne = (((NLO_PHI)?h_nlo_phi_scale_m->GetBinContent(bin_lo):h_lo_phi_scale_m->GetBinContent(bin_lo)) / ((NLO_QCD)?h_nlo_qcd_scale_m->GetBinContent(bin_lo):h_lo_qcd_scale_m->GetBinContent(bin_lo)) > 1.0);
  	  double sigma_qcd_m = 0.0;
  	  double sigma_phi_m = 0.0;
  	  double sigma_qcd_h = 0.0;
  	  double sigma_phi_h = 0.0;
  	  double sigma_qcd_l = 0.0;
  	  double sigma_phi_l = 0.0;
  	  double R_m = 0.0;
  	  double N, Nhi, Nlo;
  	  int ibin;
  	  int bin_x = bin_lo;

  	  for (ibin = bin_lo; ibin<bin_hi; ++ibin)
  	    {
  	      double binw      = (NLO_QCD)?h_nlo_qcd_scale_m->GetBinWidth(ibin):h_lo_qcd_scale_m->GetBinWidth(ibin);
  	      double val_qcd_m = (NLO_QCD)?h_nlo_qcd_scale_m->GetBinContent(ibin):h_lo_qcd_scale_m->GetBinContent(ibin);
  	      double val_phi_m = (NLO_PHI)?h_nlo_phi_scale_m->GetBinContent(ibin):h_lo_phi_scale_m->GetBinContent(ibin);
  	      double val_qcd_h = (NLO_QCD)?h_nlo_qcd_scale_h->GetBinContent(ibin):h_lo_qcd_scale_h->GetBinContent(ibin);
  	      double val_phi_h = (NLO_PHI)?h_nlo_phi_scale_h->GetBinContent(ibin):h_lo_phi_scale_h->GetBinContent(ibin);
  	      double val_qcd_l = (NLO_QCD)?h_nlo_qcd_scale_l->GetBinContent(ibin):h_lo_qcd_scale_l->GetBinContent(ibin);
  	      double val_phi_l = (NLO_PHI)?h_nlo_phi_scale_l->GetBinContent(ibin):h_lo_phi_scale_l->GetBinContent(ibin);
  	      if ( !subIntervals ||  ((val_phi_m / val_qcd_m > 1.0) == ratioGtOne) )
  		{
  		  sigma_phi_m += binw*val_phi_m;
  		  sigma_qcd_m += binw*val_qcd_m;
  		  sigma_phi_h += binw*val_phi_h;
  		  sigma_qcd_h += binw*val_qcd_h;
  		  sigma_phi_l += binw*val_phi_l;
  		  sigma_qcd_l += binw*val_qcd_l;
  		}
  	      else
  		{
// rrdtool graph cputemp.gif --start -24h --title "CPU-Temperatur" --vertical-label "Grad Celsius" DEF:cputemperatur=cputemp.rrd:temp:AVERAGE LINE1:cputemperatur#ff0000:"CPU-Temperatur"

  		  cout << endl;
  		  // cout << " x   range : " << setw(5) << h_lo_qcd_scale_m->GetBinLowEdge(bin_x);
  		  // cout << setw(5) << h_lo_qcd_scale_m->GetBinLowEdge(ibin) << endl;
  		  cout << " bin range : " << setw(5) << bin_x << setw(5) << ibin << endl;
		  // output expectation value (QCD only) in selected bin range with scale variation
		  cout << setw(25)  << "<O>[QCD]    :" << setw(15) << ( sigma_qcd_m );
		  cout << " mu_h: " << setw(15) << sigma_qcd_h - sigma_qcd_m;
		  cout << " mu_l: " << setw(15) << sigma_qcd_l - sigma_qcd_m << endl;
		  // output expectation value (QCD + PHI) in selected bin range with scale variation	  
		  cout << setw(25)  << "<O>[QCD+PHI]:" << setw(15) << ( sigma_phi_m );
		  cout << " mu_h: " << setw(15) << sigma_phi_h - sigma_phi_m;
		  cout << " mu_l: " << setw(15) << sigma_phi_l - sigma_phi_m << endl;
		  // output ratio in selected bin range with scale variation	  
		  cout << setw(25)  << " R          :" << setw(15) << ( R_m = sigma_phi_m / sigma_qcd_m );
		  cout << " mu_h: " << setw(15) << sigma_phi_h / sigma_qcd_h - R_m;
		  cout << " mu_l: " << setw(15) << sigma_phi_l / sigma_qcd_l - R_m << endl;
  		  // events due to resonant production in this Mtt range
  		  N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_ll*eps_ll;
  		  Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_ll*eps_ll;
  		  Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_ll*eps_ll;
  		  cout << setw(25) << " delta N NLO PHI (ll): "   << setw(15) << N;
  		  cout << " mu_h: " << setw(15) << Nhi - N;
  		  cout << " mu_l: " << setw(15) << Nlo - N << endl;
  		  N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_lj*eps_lj;
  		  Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_lj*eps_lj;
  		  Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_lj*eps_lj;
  		  cout << setw(25) << " delta N NLO PHI (lj): "   << setw(15) << N;
  		  cout << " mu_h: " << setw(15) << Nhi - N;
  		  cout << " mu_l: " << setw(15) << Nlo - N << endl;

  		  ratioGtOne = !ratioGtOne;
  		  sigma_phi_m = binw*val_phi_m;
  		  sigma_qcd_m = binw*val_qcd_m;
  		  sigma_phi_h = binw*val_phi_h;
  		  sigma_qcd_h = binw*val_qcd_h;
  		  sigma_phi_l = binw*val_phi_l;
  		  sigma_qcd_l = binw*val_qcd_l;		  
  		  bin_x = ibin;
  		}
  	    }
  	  cout << endl;
  	  // cout << " x   range : " << setw(5) << h_lo_qcd_scale_m->GetBinLowEdge(bin_x);
  	  // cout << setw(5) << h_lo_qcd_scale_m->GetBinLowEdge(ibin) << endl;
  	  cout << " bin range : " << setw(5) << bin_x << setw(5) << ibin << endl;
	  // output expectation value (QCD only) in selected bin range with scale variation
  	  cout << setw(25)  << "<O>[QCD]    :" << setw(15) << ( sigma_qcd_m );
  	  cout << " mu_h: " << setw(15) << sigma_qcd_h - sigma_qcd_m;
  	  cout << " mu_l: " << setw(15) << sigma_qcd_l - sigma_qcd_m << endl;
	  // output expectation value (QCD + PHI) in selected bin range with scale variation	  
  	  cout << setw(25)  << "<O>[QCD+PHI]:" << setw(15) << ( sigma_phi_m );
  	  cout << " mu_h: " << setw(15) << sigma_phi_h - sigma_phi_m;
  	  cout << " mu_l: " << setw(15) << sigma_phi_l - sigma_phi_m << endl;
	  // output ratio in selected bin range with scale variation	  
  	  cout << setw(25)  << " R          :" << setw(15) << ( R_m = sigma_phi_m / sigma_qcd_m );
  	  cout << " mu_h: " << setw(15) << sigma_phi_h / sigma_qcd_h - R_m;
  	  cout << " mu_l: " << setw(15) << sigma_phi_l / sigma_qcd_l - R_m << endl;
  	  // events due to resonant production in this Mtt range
  	  N   = (sigma_phi_m)*Lumi*Br_ll*eps_ll;//-sigma_qcd_m
  	  Nhi = (sigma_phi_h)*Lumi*Br_ll*eps_ll;//-sigma_qcd_h
  	  Nlo = (sigma_phi_l)*Lumi*Br_ll*eps_ll;//-sigma_qcd_l
  	  cout << setw(25) << " delta N NLO PHI (ll): "   << setw(15) << N;
  	  cout << " mu_h: " << setw(15) << Nhi - N;
  	  cout << " mu_l: " << setw(15) << Nlo - N << endl;
  	  N   = (sigma_phi_m)*Lumi*Br_lj*eps_lj;
  	  Nhi = (sigma_phi_h)*Lumi*Br_lj*eps_lj;
  	  Nlo = (sigma_phi_l)*Lumi*Br_lj*eps_lj;
  	  cout << setw(25) << " delta N NLO PHI (lj): "   << setw(15) << N;
  	  cout << " mu_h: " << setw(15) << Nhi - N;
  	  cout << " mu_l: " << setw(15) << Nlo - N << endl;
  	}
      else
  	{
  	  cout << "\n Range not valid, enter: " << xmin;
  	  cout << " < xlo < xhi < " << xmax << endl;
  	  cout << "\n You entered: xlo = " << xlo << "  , xhi = " << xhi << endl;
  	  break;
  	}
    }
  
  // produce the LO Phi+QCD / QCD ratios
  // if NLO PHI is available:
  // the histograms h_lo_phi_scale_x do not yet contain the QCD contribution

  // h_lo_phi_scale_m->Scale(2);
  // h_lo_phi_scale_h->Scale(2);
  // h_lo_phi_scale_l->Scale(2);
  
  if (NLO_PHI) h_lo_phi_scale_m->Add(h_lo_qcd_scale_m,1);
  h_lo_phi_scale_m->SetLineWidth(1);
  h_lo_phi_scale_m->SetLineColor(1);
  h_lo_phi_scale_m->SetLineStyle(3);

  if (NLO_PHI) h_lo_phi_scale_h->Add(h_lo_qcd_scale_h,1);
  h_lo_phi_scale_h->SetLineWidth(1);
  h_lo_phi_scale_h->SetLineColor(1);
  h_lo_phi_scale_h->SetLineStyle(3);

  if (NLO_PHI) h_lo_phi_scale_l->Add(h_lo_qcd_scale_l,1);
  h_lo_phi_scale_l->SetLineWidth(1);
  h_lo_phi_scale_l->SetLineColor(1);
  h_lo_phi_scale_l->SetLineStyle(3);

  
  // the NLO Phi+QCD / QCD ratios for mean, high, low scale
  if (NLO_PHI)
    {
      h_nlo_phi_scale_m->Divide(NLO_QCD?h_nlo_qcd_scale_m:h_lo_qcd_scale_m);  
      h_nlo_phi_scale_h->Divide(NLO_QCD?h_nlo_qcd_scale_h:h_lo_qcd_scale_h);
      h_nlo_phi_scale_l->Divide(NLO_QCD?h_nlo_qcd_scale_l:h_lo_qcd_scale_l);
    }
  else
    {
      h_lo_phi_scale_m->Divide(NLO_QCD?h_nlo_qcd_scale_m:h_lo_qcd_scale_m);  
      h_lo_phi_scale_h->Divide(NLO_QCD?h_nlo_qcd_scale_h:h_lo_qcd_scale_h);
      h_lo_phi_scale_l->Divide(NLO_QCD?h_nlo_qcd_scale_l:h_lo_qcd_scale_l);
    }

  
  // produce a scale band from the ratios above
  TH1D* scale = NLO_PHI?makeErrorHist(h_nlo_phi_scale_l,h_nlo_phi_scale_h):makeErrorHist(h_lo_phi_scale_l,h_lo_phi_scale_h);

  // retrieve canvas for drawing from root file  
  stringstream canvasName;   canvasName   << "Canvas_" << obj_num;
  TCanvas* c_scale_m = nullptr;
  file_scale_m->GetObject(canvasName.str().c_str(),c_scale_m);

  if (!c_scale_m)
    {
      cout << endl << " Error reading objects from .root files. Could not get canvas. " << endl;
      exit(0);
    }
  c_scale_m->ls();

  // retrieve lower pad from root file to draw the ratio
  stringstream upperPadName;
  upperPadName   << "upper_pad_" << obj_num;
  TPad* c_scale_m_upper_pad = (TPad*) (c_scale_m->FindObject(upperPadName.str().c_str()));

  stringstream lowerPadName;
  lowerPadName   << "lower_pad_" << obj_num;
  TPad* c_scale_m_lower_pad = (TPad*) (c_scale_m->FindObject(lowerPadName.str().c_str()));
  if (!c_scale_m_lower_pad)
    {
      cout << endl << " Error reading objects from .root files. Could not get lower pad in canvas. " << endl;
      exit(0);
    }



  
  TApplication TheApp("MyApp",&argc, argv);
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleAlign(0);

  // c_scale_m_lower_pad->ls("-m");
  // c_scale_m_lower_pad->RecursiveRemove(c_scale_m_lower_pad->GetPrimitive((NLO_PHI)?h_nlo_phi_scale_m->GetName():h_lo_phi_scale_m->GetName()));
  c_scale_m_lower_pad->cd();

  // c_scale_lower_pad->DrawFrame();
  
  // TCanvas* c = new TCanvas("canvas","canvas",0,0,1000,1000);
  // c->cd(1);

  // TLegend* leg = new TLegend(
  // 			     0.7428,
  // 			     0.8461,
  // 			     0.9559,
  // 			     0.9550
  // 			     );
  
  // scale->GetXaxis()->SetLabelSize(0.025);
  // scale->GetYaxis()->SetLabelSize(0.025);
  // scale->GetYaxis()->SetTitleSize(0.025);
  // scale->GetXaxis()->SetTitleSize(0.025);
  // scale->GetYaxis()->SetTitleOffset(1.3);
  // scale->GetXaxis()->SetTitleOffset(1.3);
  // scale->GetYaxis()->SetTitle("#Delta B  [pb/GeV]");
      
  scale->DrawCopy("same e2");
  if (NLO_PHI)
    {
      h_nlo_phi_scale_m->DrawCopy("same hist ][");
    }
  else
    {
      h_lo_phi_scale_m->SetAxisRange(340,750);
      h_lo_phi_scale_m->SetLineColor(2);
      h_lo_phi_scale_m->SetLineStyle(1);
      h_lo_phi_scale_m->SetLineWidth(2);
      h_lo_phi_scale_m->DrawCopy("same hist ][");
      // leg->AddEntry(h_lo_qcd_scale_m,"LO QCD","L");
      // leg->AddEntry(h_lo_phi_scale_m,"LO QCD + LO #Phi_{2},#Phi_{3}","L");
    }

  
  // h_lo_qcd_scale_m->Divide(h_lo_qcd_scale_m);
  // h_lo_qcd_scale_m->DrawCopy("same hist ][");
  // leg->Draw("same");
  // c->Draw();

  
  c_scale_m->ls();
  c_scale_m->Draw();


  TheApp.Run();

  return 1;
}














  // // analyze Mtt distribution below and above resonances
  // if (obj_num == 0)
  //   {
  //     int nbins = h_nlo_qcd_scale_m->GetNbinsX();
  //     int ibin  = 1;
  //     double sigma_qcd_m = 0.0;
  //     double sigma_phi_m = 0.0;
  //     double sigma_qcd_h = 0.0;
  //     double sigma_phi_h = 0.0;
  //     double sigma_qcd_l = 0.0;
  //     double sigma_phi_l = 0.0;
  //     double R_m = 0.0;
  //     double N, Nhi, Nlo;
	  
  //     // sum as long as ratio is greater/less than 1
  //     // if this changes return sigma, bin
  //     bool belowRes = true;
  //     // sum total cross section from Mtt distribution until the bin where the ratio is < 1
  //     while (ibin<=nbins)
  // 	{
  // 	  double binw      = h_nlo_qcd_scale_m->GetBinWidth(ibin);
  // 	  double val_qcd_m = h_nlo_qcd_scale_m->GetBinContent(ibin);
  // 	  double val_phi_m = h_nlo_phi_scale_m->GetBinContent(ibin);
  // 	  double val_qcd_h = h_nlo_qcd_scale_h->GetBinContent(ibin);
  // 	  double val_phi_h = h_nlo_phi_scale_h->GetBinContent(ibin);
  // 	  double val_qcd_l = h_nlo_qcd_scale_l->GetBinContent(ibin);
  // 	  double val_phi_l = h_nlo_phi_scale_l->GetBinContent(ibin);
  // 	  double bin_low   = h_nlo_phi_scale_m->GetBinLowEdge(ibin);
	  
	  
  // 	  if (val_phi_m / val_qcd_m < 1.0 && belowRes)
  // 	    {
  // 	      cout << endl << " last bin: " << setw(5) << ibin << endl;
  // 	      cout << setw(25)  << " R:" << setw(15) << ( R_m = sigma_phi_m / sigma_qcd_m );
  // 	      cout << " mu_h: " << setw(15) << sigma_phi_h / sigma_qcd_h - R_m;
  // 	      cout << " mu_l: " << setw(15) << sigma_phi_l / sigma_qcd_l - R_m << endl;
  // 	      // events due to resonant production in this Mtt range
  // 	      N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_ll*eps_ll;
  // 	      Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_ll*eps_ll;
  // 	      Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_ll*eps_ll;
  // 	      cout << setw(25) << " delta N NLO PHI (ll): "   << setw(15) << N;
  // 	      cout << " mu_h: " << setw(15) << Nhi - N;
  // 	      cout << " mu_l: " << setw(15) << Nlo - N << endl;
  // 	      N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_lj*eps_lj;
  // 	      Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_lj*eps_lj;
  // 	      Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_lj*eps_lj;
  // 	      cout << setw(25) << " delta N NLO PHI (lj): "   << setw(15) << N;
  // 	      cout << " mu_h: " << setw(15) << Nhi - N;
  // 	      cout << " mu_l: " << setw(15) << Nlo - N << endl;
      
  // 	      belowRes = false;
  // 	      sigma_phi_m = 0;
  // 	      sigma_qcd_m = 0;
  // 	      sigma_phi_h = 0;
  // 	      sigma_qcd_h = 0;
  // 	      sigma_phi_l = 0;
  // 	      sigma_qcd_l = 0;
  // 	    }
  // 	  // add bin contents to cross section in Mtt range
  // 	  if (bin_low>=360.0 && bin_low<800.0)// && fabs(val_phi / val_qcd - 1.0) > 0.005)
  // 	    {
  // 	      sigma_phi_m += binw*val_phi_m;
  // 	      sigma_qcd_m += binw*val_qcd_m;
  // 	      sigma_phi_h += binw*val_phi_h;
  // 	      sigma_qcd_h += binw*val_qcd_h;
  // 	      sigma_phi_l += binw*val_phi_l;
  // 	      sigma_qcd_l += binw*val_qcd_l;
  // 	    }
  // 	  else
  // 	    {
  // 	      cout << endl << " omitting bin " << ibin << endl;
  // 	    }
  // 	  ++ibin;
  // 	}
  //     cout << endl << " last bin: " << setw(5) << ibin << endl;
  //     cout << setw(25)  << " R:" << setw(15) << ( R_m = sigma_phi_m / sigma_qcd_m );
  //     cout << " mu_h: " << setw(15) << sigma_phi_h / sigma_qcd_h - R_m;
  //     cout << " mu_l: " << setw(15) << sigma_phi_l / sigma_qcd_l - R_m << endl;
  //     // events due to resonant production in this Mtt range
  //     N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_ll*eps_ll;
  //     Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_ll*eps_ll;
  //     Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_ll*eps_ll;
  //     cout << setw(25) << " delta N NLO PHI (ll): "   << setw(15) << N;
  //     cout << " mu_h: " << setw(15) << Nhi - N;
  //     cout << " mu_l: " << setw(15) << Nlo - N << endl;
  //     N   = (sigma_phi_m-sigma_qcd_m)*Lumi*Br_lj*eps_lj;
  //     Nhi = (sigma_phi_h-sigma_qcd_h)*Lumi*Br_lj*eps_lj;
  //     Nlo = (sigma_phi_l-sigma_qcd_l)*Lumi*Br_lj*eps_lj;
  //     cout << setw(25) << " delta N NLO PHI (lj): "   << setw(15) << N;
  //     cout << " mu_h: " << setw(15) << Nhi - N;
  //     cout << " mu_l: " << setw(15) << Nlo - N << endl;
  //   }

  
  // // analyze Delta Y distribution for assymetries
  // if (obj_num == 6)
  //   {
  //     double Np_m = 0.0;
  //     double Nm_m = 0.0;
  //     double Np_h = 0.0;
  //     double Nm_h = 0.0;
  //     double Np_l = 0.0;
  //     double Nm_l = 0.0;
  //     for (int i=1;i<=h_nlo_phi_scale_m->GetNbinsX();++i)
  // 	{
  // 	  double val_m  = h_nlo_phi_scale_m->GetBinContent(i);
  // 	  double val_h  = h_nlo_phi_scale_h->GetBinContent(i);
  // 	  double val_l  = h_nlo_phi_scale_l->GetBinContent(i);
	  
  // 	  double bin_wd = h_nlo_phi_scale_m->GetBinWidth(i);
  // 	  double bin_lo = h_nlo_phi_scale_m->GetBinLowEdge(i);
  // 	  double bin_hi = bin_lo + bin_wd;
  // 	  if (bin_lo < 0.0 && bin_hi <=0.0)
  // 	    {
  // 	      Nm_m += val_m;
  // 	      Nm_h += val_h;
  // 	      Nm_l += val_l;
  // 	    }
  // 	  else if (bin_lo >= 0.0 && bin_hi > 0.0)
  // 	    {
  // 	      Np_m += val_m;
  // 	      Np_h += val_h;
  // 	      Np_l += val_l;
  // 	    }
  // 	  else
  // 	    {
  // 	      Nm_m += val_m*(-bin_lo)/bin_wd;
  // 	      Np_m += val_m*(bin_hi)/bin_wd;
  // 	      Nm_h += val_h*(-bin_lo)/bin_wd;
  // 	      Np_h += val_h*(bin_hi)/bin_wd;
  // 	      Nm_l += val_l*(-bin_lo)/bin_wd;
  // 	      Np_l += val_l*(bin_hi)/bin_wd;
  // 	    }
  // 	}
  //     double Nc_m = 0.0;
  //     cout << endl << " Nc = " << setw(15) << (Nc_m = (Np_m-Nm_m)/(Np_m+Nm_m));
  //     cout << " " << setw(10) << (Np_h-Nm_h)/(Np_h+Nm_h);
  //     cout << " " << setw(10) << (Np_l-Nm_l)/(Np_l+Nm_l);
  //     cout << endl;
  //   }
