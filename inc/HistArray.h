
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
#include "Observables.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
\typedef Alias for the boost path class.
*/
typedef boost::filesystem::path boost_path;

/*!
  Better readable labels for the histograms in the HistArray. In the default setup there are 6 histograms, one for each of the following contribution: LO QCD, NLO QCD, LO (PHIxPHI + PHIxQCD), virtual corrections (PHIxPHI + PHIxQCD), integrated dipoles (PHIxPHI + PHIxQCD)  and real corrections (PHIxPHI + PHIxQCD).
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
\typedef This is used for input file operations: mapping of table column number to HistArray index, i.e. each entry in the map specifies which column in the file is read into which histogram in the array.
 */
typedef std::map<int,H_Index> H_IndexMap;

/*!
  Default mapping from data table column to histogram indices used for input file operations.
 */
extern const H_IndexMap imap_default;
      


// stores a histogram for Born(QCD), Born(INT+PHI), (empty), virtual, int. dipoles and real corrections
#define NHIST 6
/*!
  \brief This class serves as storage object for distributions.
  
  This class stores an array of ROOT TH1D's. In the default setup there are 6 histograms, one for each of the following contribution: LO QCD, NLO QCD, LO (PHIxPHI + PHIxQCD), virtual corrections (PHIxPHI + PHIxQCD), integrated dipoles (PHIxPHI + PHIxQCD)  and real corrections (PHIxPHI + PHIxQCD). Use one HistArray for each observable. Historgams in the array can be activated seperately, so that calls to FillOne/FillAll only affect the activated histograms. The PS objects are in charge for filling the histograms since the value of an observable depends on the respective phase space.
  \sa PhaseSpace.h
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
  //! Unique number for the HistArray instance. Derived from the static d_ID.
  int d_id;
  //! Buffer for weights from real corrections/unintegrated dipoles. Here we have the problem that large numerical values from real corrections and dipoles get spread over the histogram which makes the bin contents instable. 'NumLimit' defines what we consider a large number here, these numbers get collected in the buffer corresponding to each bin until the buffer value drops below the limit.
  double *d_buffer;

  
  //! HistArray ID counter.
  static int d_ID;
  
 public:
  //! Limit that seperates 'large numbers' from 'normal numbers'.
  static double NumLimit;
  static double NumMax;
  
  HistArray(int nbinsx,
	    double xlow,
	    double xup,
	    int mass_dim,
	    std::string const& name="",
	    std::string const& lab_x="",
	    std::string const& lab_y="",
	    bool SUMW2 = false);
  ~HistArray();

  //! Access a specific histogram in the array. Returns pointer to the respective TH1D instance.
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
  
  //! Deactivate all histograms temporarely, used for VEGAS warmup run.
  void Pause() { if (d_active != 0) std::swap(d_active_t,d_active); }
  bool IsPaused() { return !d_active; } 
  
  //! Reset flags to previous state after Pause() was called.
  void Resume() { if (d_active == 0) std::swap(d_active_t,d_active); }
  
  //! Set the description of the distribution.
  void SetName(std::string const& name) { d_name = name; }
  
  //! Get the description of the distribution, const char* for ROOT classes/functions.
  const char* GetName() const { return d_name.c_str(); }

  //! Get x-axis label, const char* for ROOT classes/functions.
  const char* GetLabX() const { return d_label_x.c_str(); }

  //! Get y-axis label, const char* for ROOT classes/functions.
  const char* GetLabY() const { return d_label_y.c_str(); }
  
  //! Get mass dimension of the associated observable.
  int GetMassDim() { return d_mass_dim; }
  
  //! Fill weight into all active histograms. Calls the respective TH1D member function on active histograms.
  /*!
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillAll(double const& x, double const& wgt) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Fill(x,wgt);} }
  
  //! Fill weight into a single histogram if activated. Calls the respective TH1D member function if the histogram is activated.
  /*!
    \param i Histogram index.
    \param x value of the observable -> specifies the bin that will be filled
    \param wgt weight to be added to the respective bin
  */
  void FillOne(H_Index i, double const& x, double const& wgt);
  void FlushBuffer(H_Index i);
  
  //! Create canvas with ratio plot and draw the histgrams.
  /*!
    \param writeToRootFile Write canvas and histograms to the ROOT file which is currently opened.
  */
  void DrawCanvas(bool writeToRootFile);
  
  
  //! Draw all histograms into the currently selected canvas. Calls the respective TH1D member function.
  /*!
    \param opt ROOT drawing options
  */
  void Draw(const char* opt="") { for (unsigned i=0;i<NHIST;i++) {d_histograms[i].Draw(opt);} }

  //! Rescale all active histograms by a factor of 1/c. Calls the respective TH1D member function.
  void Scale(double c) { for (unsigned i=0;i<NHIST;i++) { if(IsActive(i)) d_histograms[i].Scale(1.0/c); } }
  
  //! normalize all histograms such that sum(b_i) = sigma, where b_i are the bins contents and sigma is the total cross section. This function assumes equal bin widths!!!
  /*!
    \param mScale if dimensionful observables have been normalized to mass scale mScale this can be used to restore the original mass scale (provided d_mass_dim was set up correctly)
  */
  void Normalize(const double& mScale=1.0, int verb = 0);

  //! Reset all histograms. Calls the respective TH1D member function.
  /*!
    \param opt ROOT options
  */
  void Reset(const char* opt="") { for (unsigned i=0;i<NHIST;i++) {d_histograms[i].Reset(opt);} }
  
  //! Read data in column 'colNum' of table 'tabNum' in the file specified by 'path' into histogram 'i'. Previous data will be erased. The method checks if the histogram binning matches the input data. The expected file format is:
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
  //! Lines starting with '#' will be ignored as comments. Tables need to be separated with one comment line!
  /*!
    \param i Histogram index, data destination.
    \param x The path of the input file, data source.
    \param tabNum Pick table  'tabNum' from the input data in the file and copy to histogram 'i'.
    \param colNum Pick column 'colNum' from this table.
  */
  int ReadFile(H_Index i,
	       const boost_path& path,
	       int tabNum,
	       int colNum);


  //! Read data in the next table from ifstream 'file' into histogram 'i'. The IndexMap specifies which table column is read into which histogram in the array. Previous data will be erased. The method checks if the histogram binning matches the input data. The expected file format is:
  /*!
    # table 1
    ...
    lower_bin_edge   data_col_1  ...  data_col_n
    ...
    # table 2
    ...
    lower_bin_edge   data_col_1  ...  data_col_n
    ...
    ...
  */
  //! Lines starting with '#' will be ignored as comments. Tables need to be separated with one comment line!
  /*!
    \param file Input file, data source.
    \param norm Normalization factor applied to all data elements read from the ifstream.
    \param imap IndexMap, specifies which source column is copied to which histogram in the array.
    \param discardCurrentTable Jump to the next table if previous read operation has left the ifstream in the middle of a table.
  */  
  int ReadTable(
		std::ifstream& file,
		const double& norm = 1.0,
		const H_IndexMap& imap = imap_default,
		bool discardCurrentTable=false);
   
  //! Print content of the HistArray instance to the specified std::ostream.
  void Print(std::ostream& ost, int verb = 0);

  //! Shows which histogram in the array is currently active.
  void Status(std::ostream& ost);

  //! Get the unique ID of this HistArray instance.
  int GetID() const { return d_id; }
  
  // layout settings for the generated TH1D's
  static double LabelSize;
  static double TitleOffset;
  static double CanvasLeftMargin;
  static double CanvasBottomMargin;
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
//! top-polar angle distributions (Collins-Soper angle)
extern HistArray T1_Histograms;
//! antitop-polar angle distributions (Collins-Soper angle)
extern HistArray T2_Histograms;

// angular distributions (spin dependent)
//! Mtt distribution of lepton opening angle correlation
extern HistArray Dopen_Histograms;
//! Mtt distribution of CP-odd triple correlation
extern HistArray OCP1_Histograms;
//! Mtt distribution of top transverse polarization
extern HistArray B1_Histograms;
//! Mtt distribution of antitop transverse polarization
extern HistArray B2_Histograms;
//! Mtt distribution of helicity angle correlation
extern HistArray Chel_Histograms;




///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/*!
  \brief Class for differential distributions.
  
  The distribution class holds a pointer to a storage object (HistArray) and an associated observable function (OBSFnc, see Observables.h). The phase space classes serve as the 'source' of data for the distributions. The distribution class can thus be seen as the interface between these two.
  \sa  Observables.h, PhaseSpace.h
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
  double operator()(const PS_2* ps_lab, const PS_2* ps_tt) { return d_fnc(ps_lab,ps_tt);  }
  double Obs(const PS_2* ps_lab, const PS_2* ps_tt)        { return d_fnc(ps_lab,ps_tt);  }
  virtual double Avg(const PS_2*, const PS_2*)             { return 1.0; }
};


/*!
  \brief Class for differential distributions of mean values.
  
  The mean distribution is an extension of the distribution class. It stores a pointer to a second observable function (OBSFnc). Not only the bare event weight is filled into a histogram, but it is multiplied by a second observable, thereby computing the mean value of the latter as a distribution in the first observable.
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
  double Avg(const PS_2* ps_lab, const PS_2* ps_tt)        { return d_fnc_mean(ps_lab,ps_tt); }
};
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


/*!
\typedef A vector with pointers to distribution objects. This is actually the type of object the phase spaces are operating on when filling event data into distributions.
  \sa PhaseSpace.h
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
