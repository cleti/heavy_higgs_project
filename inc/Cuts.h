
/*! \file
  \brief Phase space cuts and collections of cuts. 
*/

#ifndef CUTS_H
#define CUTS_H

#include "Observables.h"
#include "PhaseSpace.h"
#include <vector>

/*!
  \brief Represents a single phase space cut.
  
  On checking the cut via operator '()', the observable 'd_obs' is evaluated and the result compared to the ranges given via 'd_intervals'. Each pair of values in the vector defines an interval [d_in[i],d_int[i+1]] for which the cut returns true if the value of 'd_obs' lies within. If an uneven number of values is given in 'd_intervals', the last value defines an interval [d_int[last], + infinity] for which the cut is accepted.
  \sa Observables.h
*/
class Cut {
 private:
  OBSFnc               d_obs;
  std::vector<double>  d_intervals;
  
  bool check_intervals();
  
 public:
  Cut(OBSFnc obs, const std::initializer_list<double>& il);
  ~Cut() {}

  bool set_intervals(const std::initializer_list<double>& il);
  bool operator()(const PS_2* ps_lab, const PS_2* ps_tt) const;
};


/*!
  \brief Collection of phase space cuts.
*/
typedef std::vector< std::shared_ptr<Cut> > CutVec; 

  /*!
    Check the whole collection of cuts given by the pointer 'cuts'. Checks for nullptr and returns true in this case. If a single cut returns false, the result is false and the remaining cuts will not be evaluated.
    \param cuts Vector of phase space cuts
    \param ps_lab Pointer to lab-frame phase space (used to evaluate 'd_obs' in the cut)
    \param ps_tt Pointer to tt-z.m.f. phase space (used to evaluate 'd_obs' in the cut)
  */
bool EvalCuts(const CutVec* cuts, const PS_2* ps_lab, const PS_2* ps_tt);



#endif
