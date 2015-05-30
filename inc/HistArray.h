
/*! \file
  \brief Histogram class for the computation of differential distributions.
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

enum H_TYPE
  {
    H_LO_QCD    =  0,
    H_NLO_QCD   =  1,
    H_LO_PHI    =  2,
    H_NLO_PHI_V =  3,
    H_NLO_PHI_ID=  4,
    H_NLO_PHI_R =  5
  };

#define F_H_ALL BOOST_BINARY(11 111 1)


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
  //! distribution name
  std::string d_name;
  //! latex label x-axis
  std::string d_label_x;
  //! latex label y-axis
  std::string d_label_y;
  //! mass dimension of the observable, needed for proper normalization, for example [M_tt]=1, [Y_t]=0
  int d_mass_dim;

  //! running histogram id
  static int d_ID;
  
 public:
  HistArray(int nbinsx,
	    double xlow,
	    double xup,
	    int mass_dim,
	    std::string const& name="",
	    std::string const& lab_x="",
	    std::string const& lab_y="",
	    bool SUMW2 = false);
  ~HistArray() {}

  //! access a specific histogram in the array
  /*!
    \param i histogram index, NO RANGE CHECK!
  */
  TH1D* operator[](unsigned i) { return &d_histograms[i]; }
  
  //! check if a specific histogram is active
  bool IsActive(unsigned i) { return d_active & (1<<i);}
  
  //! activate specific histogram
  void SetActive(unsigned i) { d_active = (1<<i); }

  //! activate all histograms
  void ActivateAll() { d_active = F_H_ALL; }

  //! deactivate all histograms
  void DeactivateAll() { d_active = 0; }
  
  //! change state of specific histogram without touching the others
  void ToggleActive(unsigned i) { d_active ^= (1<<i); }
  
  //! deactivate all histograms, used for VEGAS warmup run
  void Pause() { if (d_active != 0) std::swap(d_active_t,d_active); }
  bool IsPaused() { return !d_active; } 
  
  //! reset flags to previous state
  void Resume() { if (d_active == 0) std::swap(d_active_t,d_active); }
  
  //! set the description of the distribution
  void SetName(std::string const& name) { d_name = name; }
  
  //! get the description of the distribution, const char* for ROOT classes/functions
  const char* GetName() { return d_name.c_str(); }

  //! get mass dimension of the associated observable
  int GetMassDim() { return d_mass_dim; }
  
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
  void FillOne(H_TYPE i, double const& x, double const& wgt) { if(IsActive(i)) { d_histograms[i].Fill(x,wgt);} }
  
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
  void Normalize(const double& mScale=1.0, int verb = 0);

  
  void Print(std::ostream& ost);
  void Status(std::ostream& ost);

  // layout settings for the generated TH1D's
  static double LabelSizeY;
  static double TitleOffsetY;
  
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////

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
extern HistArray Mtt_Histograms;
//! top transverse momentum distributions
extern HistArray PT1_Histograms;
//! antitop transverse momentum distributions
extern HistArray PT2_Histograms;
//! top+antitop transverse momentum distributions
extern HistArray PT12_Histograms;
//! top rapidity distributions
extern HistArray Y1_Histograms;
//! antitop rapidity distributions
extern HistArray Y2_Histograms;
//! top-antitop rapidity difference distributions
extern HistArray DY_Histograms;

// angular distributions (spin dependent)
extern HistArray PHIT12_Histograms;
extern HistArray Dopen_Histograms;
extern HistArray OCP1_Histograms;
extern HistArray B1_Histograms;
extern HistArray Chel_Histograms;





///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
class PS_2;
typedef double (*OBSFnc)(const PS_2*);

class Distribution {
 private:
  HistArray* d_hist;
  OBSFnc     d_fnc;
    
 public:
 Distribution(HistArray* hist, OBSFnc fnc):
    d_hist(hist),
    d_fnc(fnc)
  {}

  virtual ~Distribution() {}

  HistArray*     GetHistograms()            { return d_hist; }
  double         operator()(const PS_2* ps) { return d_fnc(ps);  }
  virtual double Avg(const PS_2* ps)        { return 1.0;    }
};


class MeanDistribution: public Distribution {
 private:
  OBSFnc     d_fnc_mean;
    
 public:
 MeanDistribution(HistArray* hist, OBSFnc fnc, OBSFnc fnc_mean):
    Distribution(hist,fnc),
    d_fnc_mean(fnc_mean)
  {}

  ~MeanDistribution() {}

  double Avg(const PS_2* ps)        { return d_fnc_mean(ps); }
};
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////



#endif
