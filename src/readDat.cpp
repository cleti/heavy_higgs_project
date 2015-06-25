





#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

// ROOT header files
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>


#include "../inc/FileBrowser.h"




const int colWidth1 = 16;
const int colWidth2 = 20;


int main(int argc, char** argv)
{
  FileBrowser b("results/PG");
  ifstream file;
  while (!file.is_open())
    {
      file.open(b.browse());
    }
  std::string sbuff;
  double bin_lo;
  double lo_phi[3];
  double nlo_phi[3];
  double lo_qcd[3];
  double nlo_qcd[3];
  double tmp;
  while(1)
    {
      if (file.good())
	{
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
		{ // print comment
		  std::cout << sbuff << std::endl;
		}
	      else
		{ // data
		  double f_qcd = 1.0/20.0;
		  double f_phi = 1.0;
		  std::stringstream iss(sbuff);
		  ///////////////////////////////
		  // lower bin boundary
		  iss >> bin_lo;
		  // // LO PHI full, without int.
		  // iss >> tmp;
		  // iss >> tmp;
		  // iss >> tmp;
		  // LO QCD
		  iss >> lo_qcd[0];
		  // iss >> lo_qcd[1];
		  // iss >> lo_qcd[2];		      
		  // LO PHI full, with int.
		  iss >> lo_phi[0];
		  // iss >> lo_phi[1];
		  // iss >> lo_phi[2];
		  // // LO PHI eff., without int.
		  // iss >> tmp;
		  // iss >> tmp;
		  // iss >> tmp;
		  // // LO PHI eff., with int.
		  // iss >> tmp;
		  // iss >> tmp;
		  // iss >> tmp;
		  // // NLO PHI, without int.
		  // iss >> tmp;
		  // iss >> tmp;
		  // iss >> tmp;

		  // NLO QCD
		  iss >> nlo_qcd[0];
		  // iss >> nlo_qcd[1];
		  // iss >> nlo_qcd[2];		      
		  // NLO PHI, with int.
		  iss >> nlo_phi[0];
		  // iss >> nlo_phi[1];
		  // iss >> nlo_phi[2];


		  ///////////////////////////////
		  // print relevant data
		  cout << std::setw(colWidth2) << std::setprecision(5)  << bin_lo;
		  cout << std::setw(colWidth2) << std::setprecision(10) << f_qcd*lo_qcd[0];
		  cout << std::setw(colWidth2) << std::setprecision(10) << f_qcd*lo_qcd[0]+f_phi*lo_phi[0];
		  cout << std::setw(colWidth2) << std::setprecision(10) << f_qcd*lo_qcd[0]+f_qcd*nlo_qcd[0];
		  cout << std::setw(colWidth2) << std::setprecision(10) << f_qcd*lo_qcd[0]+f_phi*lo_phi[0]+f_qcd*nlo_qcd[0]+f_phi*nlo_phi[0] << std::endl;
		}
	    }
	}
      else
	{
	  if (file.eof())
	    {
	      std::cout << "\n End-of-File reached on input operation. " << std::endl;
	      break;
	    }
	  if (file.fail())
	    {
	      std::cout << std::endl << " Error! ";
	      std::cout << "\n Logical error on i/o operation. " << std::endl;
	      break;
	    }
	  if (file.bad())
	    {
	      std::cout << std::endl << " Error! ";
	      std::cout << "\n Read/writing error on i/o operation. " << std::endl;
	      break;
	    }
	}
    }



  return 0;
}
