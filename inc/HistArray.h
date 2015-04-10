
/*! \file
  \brief Histogram class for the computation of differential distributions
*/ 


#ifndef HISTARRAY_H
#define HISTARRAY_H

#include <iostream>
#include <vector>
#include <sstream>
#include <TH1.h>
#include <boost/utility.hpp>

#include "Global.h"
#include "Lorentz.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////


const char* cat(std::string str, int number);

#define H_LO_QCD  0
#define H_LO_PHI  1
#define H_NLO_V   3
#define H_NLO_ID  4
#define H_NLO_R   5

// stores a histogram for Born(QCD), Born(INT+PHI), (empty), virtual, int. dipoles and real corrections
#define NHIST 6
/*!
  This class stores an array of ROOT TH1D's. In the default setup there are 6 histograms, one for each contribution: LO QCD, LO PHIxPHI, LO PHIxQCD, virtual corrections, integrated dipoles and real corrections. Use one HistArray for each observable. There are some helper functions for plotting in Functions_Shared.h. The PS objects are in charge for filling the histograms since the value of an observable depends on the respective phase space.
  \sa Functions_Shared.h, PhaseSpace.h
*/
class HistArray {
 protected:
  //! the histograms
  TH1D d_histograms[NHIST];
  //! these flags indicate which histogram in the array is active, i.e. gets filled if FillAll()/FillOne() is called
  unsigned d_active;
  //! tmp. flags, used in Pause()/Resume()
  unsigned d_active_t;
  //! description of the dsitribution
  std::string d_label;
  //! mass dimension of the observable, needed for proper normalization, for example [M_tt]=1, [Y_t]=0
  int d_mass_dim;
  //! running histogram id
  static int d_ID;

 public:
  HistArray(int nbinsx,
	    double xlow,
	    double xup,
	    int mass_dim,
	    std::string const& label="",
	    bool SUMW2 = false);
  ~HistArray() {}

  //! access a specific histogram in the array
  /*!
    \param i histogram index, NO RANGE CHECK!
  */
  TH1D* operator[](unsigned i) { return &d_histograms[i]; }
  //! check if a specific histogram is active
  bool IsActive(unsigned i) { return d_active & (1<<i);}
  //! set the flags, use for example SetActive(BOOST_BINARY(000 000 1)) to activate only the first
  void SetActive(unsigned i) { d_active = i; }
  //! deactivate all histograms, used for VEGAS warmup run
  void Pause() { if (d_active != 0) std::swap(d_active_t,d_active); }
  //! reset flags to previous state
  void Resume() { if (d_active == 0) std::swap(d_active_t,d_active); }
  //! set the description of the distribution
  void SetLabel(std::string const& label) { d_label = label; }
  //! get the description of the distribution, const char* for ROOT classes/functions
  const char* GetLabel() { return d_label.c_str(); }
  //! fill weight into all active histograms
  /*!
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillAll(double const& x, double const& wgt) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Fill(x,wgt);} }
  //! fill weight into a single histogram if activated
  /*!
    \param i histogram index, NO RANGE CHECK!
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillOne(unsigned i, double const& x, double const& wgt) { if(IsActive(i)) { d_histograms[i].Fill(x,wgt);} }
  //! draw all histograms into the currently selected canvas
  /*!
    \param opt ROOT drawing options
  */
  void Draw(const char* opt="") { for (unsigned i=0;i<NHIST;i++) {d_histograms[i].Draw(opt);} }
  //! rescale all active histograms by a factor of 1/c
  void Scale(double c) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Scale(1.0/c); } }
  //! normalize all histograms such that sum(b_i) = sigma, where b_i are the bins contents and sigma is the total cross section. This function assumes equal bin widths!!!
  /*!
    \param mScale if dimensionful observables have been normalized to mass scale mScale this can be used to restore the original mass scale (provided d_mass_dim was set up correctly)
  */
  void Normalize(const double& mScale=1.0);
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\typedef alias for a vector of HistArray pointers
 */
typedef std::vector<HistArray*> HAvec;

////////////////////////////////
extern int DIST_N_bins;

extern double DIST_M_L;
extern double DIST_M_U;

extern double DIST_P_L;
extern double DIST_P_U;

extern double DIST_Y_L;
extern double DIST_Y_U;

extern double DIST_A_L;
extern double DIST_A_U;
////////////////////////////////

//! top-antitop invariant mass distributions
extern HistArray MttDistributions;
//! top transverse momentum distributions
extern HistArray PT1Distributions;
//! antitop transverse momentum distributions
extern HistArray PT2Distributions;
//! top+antitop transverse momentum distributions
extern HistArray PT12Distributions;
//! top rapidity distributions
extern HistArray Y1Distributions;
//! antitop rapidity distributions
extern HistArray Y2Distributions;
//! top-antitop rapidity difference distributions
extern HistArray DYDistributions;
//! lepton transverse momentum distributions
extern HistArray PTL1Distributions;
//! antilepton transverse momentum distributions
extern HistArray PTL2Distributions;
//! lepton-antilepton opening angle distributions
extern HistArray PHILLDistributions;
//! partonic c.m.e. distributions
extern HistArray SPartDistributions;




//! compute invariant mass of two 4-vectors k1 and k2
double obs_M12(FV const& k1, FV const& k2);
//! compute transverse momentum of a 4-vector k1
double obs_PT(FV const& k);
//! compute rapidity from a 4-vector k1
double obs_Y(FV const& k);
//! compute rapidity difference from two 4-vectors k1 and k2
double obs_DY(FV const& k1, FV const& k2);
//! compute momentum direction opening angle form two 4-vectors k1 and k2
double obs_PHI(FV const& k1, FV const& k2);



#endif
