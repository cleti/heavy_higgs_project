
/*! \file
 \brief This class provides the functionality of the GSL VEGAS integration algorithm.
*/ 

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <cstdlib> // size_t
#include <bitset>
#include <memory>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

#include "Global.h"
#include "HistArray.h"
#include "PhaseSpace.h"
#include "HiggsModel.h"

#include "LHAPDF/LHAPDF.h"


typedef unsigned long ulong;



///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*!
  This class contains the parameters to control the GSL VEGAS algorithm.
*/
class vegas_par
{
 public:
  int      verbose;
  int      calls;
  int      iterations;
  int      max_runs;
  int      num_runs;
  double   chisq_limit;
  double   chisq;
  bool     do_warmup;
  bool     grid_fixed;
  double   result;
  double   error;

  vegas_par();
  ~vegas_par();
  void Print(std::ostream& ost = std::cout);
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*!
\typedef integrand function type required by the GSL VEGAS algorithm
 */
typedef double (*GSLIFnc)(double*,size_t,void*);

/*!
  This class represents an integral to be computed by VEGAS. This comprises the integral dimension, integral limits and a pointer to the integrand function.
*/
class Integral {

 private:
  //! integral dimension
  size_t     d_Dim;
  //! lower integration limits
  double *   d_IntLimitLo;
  //! upper integration limits
  double *   d_IntLimitUp;
  //! pointer to integrand
  GSLIFnc    d_Integrand;

  //! initialize integral of dimension dim
  void Init(size_t dim);
  
 public:
  vegas_par LastRun;
    
  Integral();
  Integral(size_t dim);
  Integral(size_t dim, double IntLimitLo[], double IntLimitUp[]);
  ~Integral();

  Integral& operator=(Integral const& rhs);
  bool operator>(size_t rhs)  { return d_Dim>rhs; }
  bool operator<(size_t rhs)  { return d_Dim>rhs; }
  bool operator>=(size_t rhs) { return d_Dim>=rhs; }
  bool operator<=(size_t rhs) { return d_Dim>=rhs; }
  bool operator==(size_t rhs) { return d_Dim==rhs; }
  bool operator!=(size_t rhs) { return d_Dim!=rhs; }

  //! check if index i is covered by the current integral size
  bool InRange(size_t i) { return (i>=0 && i<d_Dim); }
  //! get lower integration limit in direction i 
  double Lo(size_t i) { if (InRange(i)) return d_IntLimitLo[i]; else return 0.0; }
  //! get upper integration limit in direction i 
  double Up(size_t i) { if (InRange(i)) return d_IntLimitUp[i]; else return 0.0; }
  //! get pointer to lower integration limits
  double* Lo() { return d_IntLimitLo; }
  //! get pointer to upper integration limits
  double* Up() { return d_IntLimitUp; }
  //! get integral dimension
  unsigned GetDim() { return d_Dim; }
  //! set integrand
  void SetIntegrand(GSLIFnc Integrand) { d_Integrand = Integrand; }
  GSLIFnc GetIntegrand() { return d_Integrand; }
  void PrintLimits(std::ostream& ost = std::cout);
  void Print(std::ostream& ost = std::cout);
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*!
  some intergration intervals that are used in my program
*/
namespace IntLimits
{
  // 2->1 with PDF and convolution for int. dipoles
  extern double Int_lo_2_1_pdf_x[2];
  extern double Int_up_2_1_pdf_x[2];
  
  // 2->2 without PDF
  extern double Int_lo_2_2[2];
  extern double Int_up_2_2[2];
  
  // 2->2 with PDF
  extern double Int_lo_2_2_pdf[3];
  extern double Int_up_2_2_pdf[3];

  // .. including convolution for int. dipoles
  extern double Int_lo_2_2_pdf_x[4];
  extern double Int_up_2_2_pdf_x[4];

  // 2->3 with PDF
  extern double Int_lo_2_3_pdf[6];
  extern double Int_up_2_3_pdf[6];

  // 2->2 -> 3 + 3 with PDF
  // including convolution for integrated dipoles
  extern double Int_lo_2_6_pdf_x[14];
  extern double Int_up_2_6_pdf_x[14];

  // 2->3 -> 1+ 3 + 3 with PDF
  // including convolution for integrated dipoles
  extern double Int_lo_2_7_pdf[16];
  extern double Int_up_2_7_pdf[16];
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*!
  This class contains a set of parameters which is handed down to the integrand called by the GSL VEGAS algorithm. In particular the model parameters contained in the HiggsModel instance are needed by the functions that evaluate the matrix elements.
*/
class integrand_par
{
 public:
  //! pointer to object with model parameters
  HiggsModel*              higgs_model;
  //! hadronic c.m.e.
  double                   s_hadr;
  //! pointer to phase space
  PS_2*                    ps;
  //! pointer to the parton distribution function
  LHAPDF::PDF*             pdf;
  //! pointer to the GSL VEGAS state
  gsl_monte_vegas_state*   v_state;
  //! pointer to histograms with distributions
  std::vector<HistArray*>* distributions;
  //! switch to enable distribution collection
  bool                     collect_dist;
  //! include the decays of top/antitop
  bool                     tDecay;
  //! output stream
  std::ostream&            ost;
  //! evaluation flags
  ulong                    eval_flags;
  //! rescaling factor |full ggH|^2 / |eff. ggH|^2 
  double                   K;
  //! delete ps and pdf on destruction
  bool                     cleanup;           
  
  integrand_par(std::ostream& os = std::cout);
  ~integrand_par();
  //! compute the weight of the VEGAS algorithm at the current point
  double cmp_v_weight();
  //! set ps pointer, the former PS instance will be deleted!
  void SetPS(PS_2* psn) { if (ps!=nullptr) delete ps; ps = psn; }
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*!
  Wrapper class for the GSL VEGAS algorithm. Does the integration, checking integration parameters, output.
*/
class Integrator {

 private:
  std::ostream&          d_Ost;
  //! GSL integrand function 
  gsl_monte_function     d_GSLIntegrand;
  //! GSL random number generator
  gsl_rng*               d_GSLRng;
  //! GSL VEGAS state
  gsl_monte_vegas_state* d_GSLState;
  //! if an external state was passed  it will not be deleted after finishing the integration
  bool d_externalGSLState;

  //! check integrand paramters and vegas parameters for consistency
  int Init(Integral& integral, integrand_par& ip, vegas_par& vp);

 public:
  Integrator(std::ostream& ost = std::cout);
  ~Integrator();

  void SetState(gsl_monte_vegas_state* extGSLState) { d_GSLState = extGSLState; d_externalGSLState = true; }
  void DropState() { d_GSLState = nullptr; d_externalGSLState = false; }
  /*!
    invoces the GSL VEGAS algorithm with parameters set in vp, parameters in ip are handed down to the integrand
  */
  void Integrate(Integral& integral, integrand_par& ip, vegas_par& vp);
  void Reset();
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////










#endif


