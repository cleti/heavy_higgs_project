

/*! \file
  \brief Some shared helper functions.
*/ 

#ifndef FUNCTIONS_SHARED_H
#define FUNCTIONS_SHARED_H


#include <bitset>

#include "Global.h"
#include "Makros.h"
#include "HistArray.h"
#include "HiggsModel.h"

#include <sys/syscall.h>


// GSL header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_dilog.h>

#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TH2.h>
#include <TFile.h>
#include <TLegend.h>

#include "LHAPDF/LHAPDF.h"



///// ARG /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// COREUTILS header files
#include <getopt.h>

extern char *program_name;
extern struct option const longopts[];


//! Structure to represent the command line options given by the user.
struct opt {
  int int_flags;
  int n_calls;
  double ren_scale;
  double ren_scale_mult;  
  double cme;
  double as;
  double mt;
  double mb;  
  double vh;   
  int tech_cut;
  int precision;
  bool dist;
  bool tdecay;
  bool useK;
  bool logfile;
  std::string logfile_path;
  bool rootfile;
  std::string rootfile_path;
  bool qcdfile;
  std::string qcdfile_path;  
  int verb_level;
};

extern struct opt g_options;

//! Parse the command line options given by the user. The results are stored in g_options.
void parse_arguments(int argc, char** argv, struct opt& options, HiggsModel& hm);
//! Print usage information.
void usage (int status);
//! Print program version
void version (int status);

int read_qcd_data(DistVec* dist, std::string const& path, std::string const& filename);
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



// no system call wrapper in glibc
inline pid_t gettid(void)
{
  return syscall(SYS_gettid);
}



////////////////////////////////////////////////////////////////////////////////////////////////
// helper functions for drawing ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


struct CanvasPtr
{
  TCanvas * c;
  TPad* p1_1;
  TPad* p1_2;
};

CanvasPtr MakeCanvas(const HistArray& histograms,
		     int width=1000,
		     int height=1000,
		     double left_marg_scale=1.0,
		     double bottom_marg_scale=1.0);

void DrawDistribution(
		      CanvasPtr& canvas,
		      HistArray& histograms,
		      const HiggsModel& THDM,
		      bool WRITE=true,
		      bool NLO=true,
		      int CD=0);

void SetRatioPlot(TH1D* hist);
void TH1RelDiff(TH1D* h1, const TH1D* h2);

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

void match_AlphaS_f5_to_f6(double const& mu2,const int VERB=1);

void set_spins_in_tt_zmf(FV const& k1,FV const& k2, FV& s1, FV& s2,FV const& s1_r,FV const& s2_r);


////////////////////////////////////////////////////////////////////////////
// complex Phi propagator denominator
////////////////////////////////////////////////////////////////////////////
double const& DenS2(double const& s, 
		    double const& m, 
		    double const& g);
inline double const& DenS2(FV& p1, 
		    double const& m, 
		    double const& g)
{
  return DenS2(MSQ(p1),m,g);
}
inline double const& DenS2(FV&& p1, 
		    double const& m, 
		    double const& g)
{
  return DenS2(MSQ(p1),m,g);
}
c_double const& DenS(double const& s, 
		     double const& m, 
		     double const& g);
inline c_double const& DenS(FV& p1, 
		     double const& m, 
		     double const& g)
{
  return DenS(MSQ(p1),m,g);
}
inline c_double const& DenS(FV&& p1, 
		     double const& m, 
		     double const& g)
{
  return DenS(MSQ(p1),m,g);
}
inline double Den(FV const& p1, double const& m)
{
  return 1.0/(MSQ(p1)-m*m);
}



c_double LN(double const& b, int C=+1);

double EPS_(FV const& k1, FV const& k2, FV const& k3, FV const& k4);

TH1D* fabs(TH1D* hist);


extern TH1D g_hist_uid_wgts;
extern TH1D g_hist_r_wgts;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////



#endif
