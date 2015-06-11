





#include "boost/filesystem.hpp" 
#include <iostream>
#include <fstream>
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






class FileBrowser
{
private:
  boost::filesystem::path               cdir_path;
  boost::filesystem::directory_iterator cdir_itr;
  boost::filesystem::directory_iterator end_itr;
  std::vector<boost::filesystem::directory_entry> cdir_content;
  std::string ext_filter;
  
  int ls_cdir();
  std::string get_entry_path(int entry_num);
  int select_entry(int entry_num);
  int apply_filter(const boost::filesystem::path &path);
  
public:
  FileBrowser(const std::string &base_path = ".");
  ~FileBrowser() {}

  std::string browse();
  int set_cdir(const boost::filesystem::path &path);
  
};

int FileBrowser::ls_cdir()
{
  cdir_content.clear();
  int cnt = 0;
  for ( cdir_itr = boost::filesystem::directory_iterator(cdir_path);
	cdir_itr != end_itr;
	++cdir_itr )
    {
      if(is_regular_file(*cdir_itr) || is_directory(*cdir_itr) )
	{
	  if (apply_filter(cdir_itr->path())) continue;
	  std::cout << std::setw(5) << cnt << std::setw(20) << cdir_itr->path().filename() << std::endl;
	  cdir_content.push_back(*cdir_itr);
	  ++cnt;
	}
    }
  
  ext_filter.clear();
  return cnt;
}


int FileBrowser::select_entry(int entry_num)
{
  cdir_itr = boost::filesystem::directory_iterator(cdir_path);
  if (entry_num<0) return 0;
  for (int i=0;i<entry_num;++i)
    {
      ++cdir_itr;
      if (cdir_itr == end_itr) return 0;
    }
  return 1;
}

int FileBrowser::apply_filter(const boost::filesystem::path &path)
{
  std::string str = path.filename().string();
  if (str.front() == '.') return 1;
  if (str.back()  == '~') return 1;
  std::string ext = path.extension().string();
  if (!ext_filter.empty())
    {
      if (!(ext_filter == ext)) return 1;
    }
  return 0;
}


FileBrowser::FileBrowser(const std::string &base_path)
{
  if (!set_cdir(base_path))
    {
      std::cout << std::endl << " FileBrowser:  Path '" << base_path << "' not found." << std::endl;
    }
}

int FileBrowser::set_cdir(const boost::filesystem::path &path)
{
  if(boost::filesystem::exists(path))
    {
      cdir_path = path;
      return 1;
    }
  else
    {
      cdir_path = boost::filesystem::current_path();
      return 0;
    }
}

std::string FileBrowser::browse()
{
  int nentries = 0;
  int entry_num = -1;
  std::string user_input;
  std::stringstream user_input_stream;
  
  boost::filesystem::directory_iterator entry;
  while(1)
    {
      user_input_stream.clear();
      user_input_stream.flush();
      std::cout << std::endl << " Content of: " << cdir_path << ":";
      if (!ext_filter.empty())
	{
	  std::cout << "  [filtering " << ext_filter << " files ]";
	}
      std::cout <<  std::endl;

      // list entries in the current directory
      nentries = ls_cdir();

      std::cout << std::endl << " Choose entry number [u:up, q:quit, .ext:set filter]: " << std::endl; 
      std::cout << " >> ";
      std::cin >> user_input;
      user_input_stream.str(user_input);

      char tc = ' ';
      // inspect the first character of user input (ignore spaces)
      while (tc==' ')
	{
	  user_input_stream.get(tc);
	}
      user_input_stream.unget();
      if (tc == '.')
	{
	  // set extension filter
	  user_input_stream >> ext_filter;
	  if (ext_filter == ".") ext_filter.clear();
	  continue;
	}
      else if (tc == 'u' || tc == 'U')
	{
	  // go to parent dir.
	  cdir_path /= "..";// cdir_path.parent_path();
	  continue;
	}
      else if (tc == 'q' || tc == 'Q')
	{
	  // quit
	  std::cout << std::endl << " FileBrowser:  Exit." << std::endl;
	  return "";
	}
      else if (tc >= '0' && tc <= '9')
	{
	  // extract entry number
	  user_input_stream >> entry_num;
	}
      else
	{
	  std::cout << std::endl << " Invalid input " << user_input << "." << std::endl;
	  continue;
	}
      
      if (entry_num>=0 && entry_num<nentries)
	{
	  // regular file selected
	  if (boost::filesystem::is_regular_file(cdir_content[entry_num]))
	    {
	      // the path to the regular file will be returned
	      break;
	    }
	  // directory selected
	  else if (boost::filesystem::is_directory(cdir_content[entry_num]))
	    {
	      if (!set_cdir(cdir_content[entry_num].path()))
		{
		  std::cout << std::endl << " Could not select new cdir: " << cdir_content[entry_num].path() << std::endl;
		}
	      continue;
	    }
	}
      else
	{
	  std::cout << std::endl << " Could not select entry " << entry_num << "." << std::endl;
	  entry_num = -1;
	}
    }
  return cdir_content[entry_num].path().string();
}




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
		  double f_qcd = 2.0;
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
