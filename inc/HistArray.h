
/*! \file
  \brief Histogram and distribution classes for the computation of differential observables.
*/ 


#ifndef HISTARRAY_H
#define HISTARRAY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <TH1.h>
#include <boost/utility.hpp>

#include "Global.h"
#include "Lorentz.h"
#include "FileBrowser.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\typedef Alias for the boost path class.
 */
typedef boost::filesystem::path boost_path;

/*!
 Used to specify histograms in the HistArray class. In the default setup there are 6 histograms, one for each contribution: LO QCD, NLO QCD, LO PHIxPHI + LO PHIxQCD, virtual corrections, integrated dipoles and real corrections.
*/
enum H_Index
  {
    H_LO_QCD    =  0,
    H_NLO_QCD   =  1,
    H_LO_PHI    =  2,
    H_NLO_PHI_V =  3,
    H_NLO_PHI_ID=  4,
    H_NLO_PHI_R =  5
  };

/*!
\typedef For input file operations: mapping of table column number to HistArray index.
 */
typedef std::map<int,H_Index> H_IndexMap;

/*!
  Default mapping from data table column to histogram indices used for input file operations.
 */
extern const H_IndexMap imap_default;
      


// stores a histogram for Born(QCD), Born(INT+PHI), (empty), virtual, int. dipoles and real corrections
#define NHIST 6
/*!
  This class stores an array of ROOT TH1D's. In the default setup there are 6 histograms, one for each contribution: LO QCD, NLO QCD, LO PHIxPHI + LO PHIxQCD, virtual corrections, integrated dipoles and real corrections. Use one HistArray for each observable. There are some helper functions for plotting in Functions_Shared.h. The PS objects are in charge for filling the histograms since the value of an observable depends on the respective phase space.
  \sa Functions_Shared.h
*/
class HistArray {
 protected:
  //! The histograms
  TH1D d_histograms[NHIST];
  //! These flags indicate which histogram in the array is active, i.e. gets filled if FillAll()/FillOne() is called
  unsigned d_active;
  //! tmp. flags, used in Pause()/Resume()
  unsigned d_active_t;
  //! Distribution name
  std::string d_name;
  //! Latex label x-axis
  std::string d_label_x;
  //! Latex label y-axis
  std::string d_label_y;
  //! Mass dimension of the observable, needed for proper normalization, for example [M_tt]=1, [Y_t]=0
  int d_mass_dim;

  //! Histogram ID
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

  //! Access a specific histogram in the array.
  /*!
    \param i Histogram index.
  */
  TH1D* operator[](H_Index i) { return &d_histograms[i]; }
  // TH1D* operator[](int    i) { return &d_histograms[i]; }
  
  //! Check if a specific histogram is active.
  bool IsActive(unsigned i) { return d_active & (1<<i);}
  
  //! Activate only one specific histogram.
  void SetActive(unsigned i) { d_active = (1<<i); }

  //! Activate all histograms.
  void ActivateAll() { d_active = BOOST_BINARY(11 111 1); }

  //! Deactivate all histograms.
  void DeactivateAll() { d_active = 0; }
  
  //! Change state of specific histogram without touching the others.
  void ToggleActive(unsigned i) { d_active ^= (1<<i); }
  
  //! Deactivate all histograms, used for VEGAS warmup run.
  void Pause() { if (d_active != 0) std::swap(d_active_t,d_active); }
  bool IsPaused() { return !d_active; } 
  
  //! Reset flags to previous state.
  void Resume() { if (d_active == 0) std::swap(d_active_t,d_active); }
  
  //! Set the description of the distribution.
  void SetName(std::string const& name) { d_name = name; }
  
  //! Get the description of the distribution, const char* for ROOT classes/functions.
  const char* GetName() { return d_name.c_str(); }

  //! Get mass dimension of the associated observable.
  int GetMassDim() { return d_mass_dim; }
  
  //! Fill weight into all active histograms.
  /*!
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillAll(double const& x, double const& wgt) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Fill(x,wgt);} }
  
  //! Fill weight into a single histogram if activated.
  /*!
    \param i Histogram index.
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillOne(H_Index i, double const& x, double const& wgt) { if(IsActive(i)) { d_histograms[i].Fill(x,wgt);} }
  
  //! Draw all histograms into the currently selected canvas.
  /*!
    \param opt ROOT drawing options
  */
  void Draw(const char* opt="") { for (unsigned i=0;i<NHIST;i++) {d_histograms[i].Draw(opt);} }

  //! Rescale all active histograms by a factor of 1/c.
  void Scale(double c) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Scale(1.0/c); } }
  
  //! normalize all histograms such that sum(b_i) = sigma, where b_i are the bins contents and sigma is the total cross section. This function assumes equal bin widths!!!
  /*!
    \param mScale if dimensionful observables have been normalized to mass scale mScale this can be used to restore the original mass scale (provided d_mass_dim was set up correctly)
  */
  void Normalize(const double& mScale=1.0, int verb = 0);

  //! Reset all histograms.
  /*!
    \param opt ROOT  options
  */
  void Reset(const char* opt="") { for (unsigned i=0;i<NHIST;i++) {d_histograms[i].Reset(opt);} }
  
  //! Read data from file into histogram 'i'. The method checks if the histogram binning matches the input data. Expecting format:
  /*!
    # table 1
    ...
    lower_bin_edge |  data_col_1 | ... | data_col_n
    ...
    # table 2
    ...
    lower_bin_edge |  data_col_1 | ... | data_col_n
    ...
    ...
  */
  //! Lines starting with '#' will be ignored. Tables need to be separated with one comment line!
  /*!
    \param i Histogram index.
    \param x The path of the input file.
    \param tabNum Pick table  'tabNum' from the input data in the file and copy to histogram 'i'.
    \param colNum Pick column 'colNum' from this table.
  */
  int ReadFile(H_Index i, const boost_path& path, int tabNum, int colNum);
  int ReadTable(
		std::ifstream& file,
		const double& norm = 1.0,
		const H_IndexMap& imap = imap_default,
		bool discardCurrentTable=false);
   
  
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
extern HistArray B2_Histograms;
extern HistArray Chel_Histograms;





///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
class PS_2;
/*!
\typedef Observable function type.
 */
typedef double (*OBSFnc)(const PS_2*);

/*!
  The distribution class holds a pointer to a storage object (HistArray) and an associated observable function (OBSFnc). Some predefined OBSFnc's can be found in PhaseSpace.h.
  \sa PhaseSpace.h
*/
class Distribution {
 private:
  //! Pointer to the HistArray.
  HistArray* d_hist;
  //! Pointer to observable function.
  OBSFnc     d_fnc;
    
 public:
 Distribution(HistArray* hist, OBSFnc fnc):
    d_hist(hist),
    d_fnc(fnc)
  {}

  virtual ~Distribution() {}

    
  //! Return pointer to the HistArray.
  HistArray*     GetHistograms()            { return d_hist; }
  const char*    GetName()                  { return d_hist->GetName(); }
  //! Evaluate the observable at the given phase space point.
  /*!
    \param ps Phase space point
  */  
  double         operator()(const PS_2* ps) { return d_fnc(ps);  }
  virtual double Avg(const PS_2* ps)        { return 1.0;    }
};


/*!
  The mean distribution is an extension of the distribution class. It stores a pointer to a second observable function (OBSFnc). The class serves to compute distributions of the mean value of this second observable.
  \sa PhaseSpace.h
*/
class MeanDistribution: public Distribution {
 private:
  //! Pointer to the observable which will be averaged.
  OBSFnc     d_fnc_mean;
    
 public:
 MeanDistribution(HistArray* hist, OBSFnc fnc, OBSFnc fnc_mean):
    Distribution(hist,fnc),
    d_fnc_mean(fnc_mean)
  {}

  ~MeanDistribution() {}

  //! Evaluate the observable at the given phase space point.
  /*!
    \param ps Phase space point
  */    
  double Avg(const PS_2* ps)        { return d_fnc_mean(ps); }
};
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


/*!
\typedef A vector with pointers to distributions.
 */
typedef std::vector< std::shared_ptr<Distribution> > DistVec; 



int ReadData(
	     std::string path,
	     DistVec::iterator it_start,
	     DistVec::const_iterator it_end,
	     const double& norm = 1.0,
	     const H_IndexMap& imap = imap_default
	     );



#endif
