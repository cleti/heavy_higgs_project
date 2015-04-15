

#include "../inc/Integrator.h"




///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
vegas_par::vegas_par():
  verbose(1),
  calls(10000000),
  iterations(5),
  max_runs(10),
  num_runs(0),
  chisq_limit(0.2),
  chisq(0.0),
  do_warmup(true),
  grid_fixed(false),
  result(0.0),
  error(0.0)
{}
vegas_par::~vegas_par()
{}
void vegas_par::Print(std::ostream& ost)
{
  ost << std::endl;
  ost << "  verbose     = " << std::setw(10) << verbose << std::endl;
  ost << "  calls       = " << std::setw(10) << calls << std::endl;
  ost << "  iterations  = " << std::setw(10) << iterations << std::endl;
  ost << "  max_runs    = " << std::setw(10) << max_runs << std::endl;
  ost << "  num_runs    = " << std::setw(10) << num_runs << std::endl;
  ost << "  chisq_limit = " << std::setw(10) << chisq_limit << std::endl;
  ost << "  chisq       = " << std::setw(10) << chisq << std::endl;
  ost << "  do_warmup   = " << std::setw(10) << do_warmup << std::endl;
  ost << "  grid_fixed  = " << std::setw(10) << grid_fixed << std::endl;
  ost << "------------------------------------------------------------------" << std::endl;
  int prec = (int)log10(result/error)+3;
  ost << "  result      = " << std::setw(10) << std::setprecision(prec) << result << std::endl;
  ost << "  error       = " << std::setw(10) << std::setprecision(3)    << error << std::endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void Integral::Init(size_t dim)
{
  if (d_IntLimitLo != nullptr)
    {
      delete []d_IntLimitLo;
      d_IntLimitLo = nullptr;
    }
  if (d_IntLimitLo != nullptr)
    {
      delete []d_IntLimitUp;
      d_IntLimitUp = nullptr;
    }
  if (dim<=0)
    {
      d_Dim = 0;
    }
  else
    {
      d_IntLimitLo = new double[dim];
      d_IntLimitUp = new double[dim];
      d_Dim = dim;
    }
}

Integral::Integral() :
d_Dim(0),
d_IntLimitLo(nullptr),
d_IntLimitUp(nullptr),
d_Integrand(nullptr)
{}

Integral::Integral(size_t dim) :
  Integral()
{
  Init(dim);
}

Integral::Integral(size_t dim, double IntLimitLo[], double IntLimitUp[]) :
  Integral(dim)
{
  for (size_t i=0;i<d_Dim;++i)
    {
      d_IntLimitLo[i] = IntLimitLo[i];
      d_IntLimitUp[i] = IntLimitUp[i];
    }
}
Integral::~Integral()
{
  Init(0);
}

Integral& Integral::operator=(Integral const& rhs)
{
  Init(rhs.d_Dim);
  for (size_t i=0;i<d_Dim;++i)
    {
      d_IntLimitLo[i] = rhs.d_IntLimitLo[i];
      d_IntLimitUp[i] = rhs.d_IntLimitUp[i];
    }
  return *this;
}

void Integral::PrintLimits(std::ostream& ost)
{
  // print integral limits
  ost << std::endl << " Integration intervals:"  << std::endl << "  [";
  for (size_t i=0;i<d_Dim;++i) {
    ost << std::setw(6) << std::setprecision(4) << d_IntLimitLo[i] << " | ";
  }
  ost << "]" << std::endl << "  [";
  for (size_t i=0;i<d_Dim;++i) {
    ost << std::setw(6) << std::setprecision(4) << d_IntLimitUp[i] << " | ";
  }
  ost << "]" << std::endl;
}

void Integral::Print(std::ostream& ost)
{
  ost << std::endl;
  ost << "==================================================================" << std::endl;
  ost << " integrand = " << (void*)d_Integrand << std::endl;
  PrintLimits();   // integration boundaries
  LastRun.Print(); // results from last run
  ost << "==================================================================" << std::endl;
  
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
integrand_par::integrand_par(std::ostream& os):
  higgs_model(nullptr),
  s_hadr(0.0),
  ps(nullptr),
  pdf(nullptr),
  v_state(nullptr),
  distributions(nullptr),
  collect_dist(false),
  tDecay(false),
  ost(os),
  eval_flags(0),
  K(1.0),
  cleanup(true)
{}
    
integrand_par::~integrand_par()
{
  if (cleanup)
    {
      // ip holds only a copy of the GSLState pointer
      // it is deleted by the Integrator or should be handled by the user if external
      if (ps != nullptr)
      	{
      	  delete ps;
      	}
      if (pdf != nullptr)
      	{
      	  delete pdf;
      	}
      if (distributions != nullptr)
      	{
      	  delete distributions;
      	}
    }
}

#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])
double integrand_par::cmp_v_weight() {
  if (v_state == nullptr) {
    ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
    ost << " ERROR: no integrand specified!" << std::endl;
    return 0.0;
  }
  // re-calculate VEGAS weight if we have a new event
  // jac = tot. integration volume / (no. of evaluations per grid-bin)
  double jac = v_state->jac; 

  // bin of VEGAS grid where the last function call lies
  int*   bin = v_state->bin; 

  // VEGAS weight = volume of grid-bin
  double   vol = 1.0;
  unsigned dim = v_state->dim;
  for (unsigned j = 0; j < dim; ++j)
    {
      double bin_width;
      unsigned k = bin[j];
      if (k == 0)
	{
	  bin_width = COORD(v_state, 1, j);
	}
      else
	{
	  bin_width = COORD(v_state, k + 1, j) - COORD(v_state, k, j);
	}
      vol *= bin_width;
    }
  return jac*vol;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// integration intervals
namespace IntLimits
{
  using namespace Constants;

  // 2->1 with PDF
  // --------------------------- x1 , x 
  double Int_lo_2_1_pdf_x[2] = { 0.0, 0.0};
  double Int_up_2_1_pdf_x[2] = { 1.0, 1.0};
  
  // 2->2 without PDF
  // --------------------- s  ,  y
  double Int_lo_2_2[2] = { 0.0, -1.0};
  double Int_up_2_2[2] = { 1.0,  1.0};
  
  // 2->2 with PDF
  // ------------------------- x1 , x2 ,  y
  double Int_lo_2_2_pdf[3] = { 0.0, 0.0, -1.0};
  double Int_up_2_2_pdf[3] = { 1.0, 1.0,  1.0};

  // including convolution for integrated dipoles
  // --------------------------- x1 , x2 ,  y  , x
  double Int_lo_2_2_pdf_x[4] = { 0.0, 0.0, -1.0, 0.0};
  double Int_up_2_2_pdf_x[4] = { 1.0, 1.0,  1.0, 1.0};

  // 2->3 with PDF
  // ------------------------ x1 , x2 ,  y  ,   M12,  y12, phi12
  double Int_lo_2_3_pdf[6] = {0.0, 0.0, -1.0,   2.0, -1.0, 0.0  };
  double Int_up_2_3_pdf[6] = {1.0, 1.0,  1.0, 100.0,  1.0, TwoPi};

  // 2->2 -> 3 + 3 with PDF
  // including convolution for integrated dipoles
  // ------------------------------------- [   2->2 ]  [          1 -> 3 (a)         ][        1->3  (b)           ] 
  // ---------------------------- x1 , x2 ,  y  , x  ,  y3 , phi3 , M12,  y12, phi12,  y3 , phi3 , M12,  y12, phi12
  double Int_lo_2_6_pdf_x[14] = { 0.0, 0.0, -1.0, 0.0, -1.0,   0.0, 0.0, -1.0,   0.0, -1.0,   0.0, 0.0, -1.0,   0.0};
  double Int_up_2_6_pdf_x[14] = { 1.0, 1.0,  1.0, 1.0,  1.0, TwoPi, 1.0,  1.0, TwoPi,  1.0, TwoPi, 1.0,  1.0, TwoPi};

  // 2->3 -> 1+ 3 + 3 with PDF
  // including convolution for integrated dipoles
  // ----------------------------------- [          2->3         ] [          1 -> 3 (a)         ][        1->3  (b)           ] 
  // -------------------------- x1 , x2 ,  y  ,  M12 ,  y12, phi12,  y3 , phi3 , M12,  y12, phi12,  y3 , phi3 , M12,  y12, phi12
  double Int_lo_2_7_pdf[16] = { 0.0, 0.0, -1.0,   2.0, -1.0, 0.0  , -1.0,   0.0, 0.0, -1.0,   0.0, -1.0,   0.0, 0.0, -1.0,   0.0};
  double Int_up_2_7_pdf[16] = { 1.0, 1.0,  1.0, 100.0,  1.0, TwoPi,  1.0, TwoPi, 1.0,  1.0, TwoPi,  1.0, TwoPi, 1.0,  1.0, TwoPi};
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
int Integrator::Init(Integral& integral, integrand_par& ip, vegas_par& vp)
{
  int verb    = vp.verbose;
  size_t dim  = integral.GetDim();

  if (!ip.eval_flags)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: no eval_flags set" << std::endl;
      return 0;
    }


  if (!ip.higgs_model)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: no model parameters found" << std::endl;
      return 0;
    }
  
  if (verb>0)
    {
      d_Ost << "\n ====================================================================== ";
      d_Ost << "\n ======================= INTEGRATION INIT ============================= ";
    }

  if (verb>-1)
    {
      d_Ost << std::endl;
      d_Ost << std::endl << " EVAL_FLAGS = " << std::bitset<32>(ip.eval_flags).to_string() << std::endl;
    }
  /////////////////////////////////////////////////////////////////////////////////////
  // check IntegralRange
  if (dim<=0)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: integral not set up properly, dim = " << dim << std::endl;
      return 0;
    }
  if (verb>0)
    {
      integral.PrintLimits();
    }
    
  // check if the integrand fnc is specified
  if (integral.GetIntegrand() == nullptr)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: no integrand function found, nullptr received!" << std::endl;
      return 0;
    }
  /////////////////////////////////////////////////////////////////////////////////////
  
  // set up GSL Rng
  if (d_GSLRng==nullptr)
    {
      gsl_rng_env_setup();
      gsl_rng_default_seed = time(0);
      d_GSLRng = gsl_rng_alloc(gsl_rng_taus);
    }
  if (verb>0)
    {
      d_Ost << "\n  rng type   : " << gsl_rng_name(d_GSLRng);
      d_Ost << "\n  rng seed   : " << gsl_rng_default_seed << std::endl;
    }
    
  // set up GSL state / initialize external state
  if (d_externalGSLState)
    {
      if (d_GSLState==nullptr)
	{
	  d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
	  d_Ost << " ERROR: external state is NULL!" << std::endl;
	  return 0;
	}
      gsl_monte_vegas_init (d_GSLState);
      // check if the external state has the correct size
      if (dim != d_GSLState->dim)
	{
	  d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
	  d_Ost << " ERROR: External GSL state dimension does not match integral dimension!" << std::endl;
	  d_Ost << " dim(integral) = " << dim << "   dim(state) = " << d_GSLState->dim << std::endl;
	  return 0;
	}
      if (verb>0) d_Ost << "\n  (using external v_state " << d_GSLState << ") \n";
    }
  else
    {
      // always allocate new state if it is not from ext
      if (d_GSLState!=nullptr)
	{
	  if (verb>0) d_Ost << "\n  (deallocating v_state " << d_GSLState << ") \n";
	  gsl_monte_vegas_free(d_GSLState);
	  d_GSLState = nullptr;
	}
      d_GSLState = gsl_monte_vegas_alloc(dim);
      if (verb>0) d_Ost << "\n  (new v_state allocated " << d_GSLState << ") \n";
    }
  // copy pointer to integrand_par struct, so that it is available in the integrand!
  ip.v_state = d_GSLState;

  // set up the GSL Integrand struct
  d_GSLIntegrand = {integral.GetIntegrand(), dim, (void*)&ip};

  if (ip.ps == nullptr)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: no ps found, nullptr received!" << std::endl;
      return 0;
    }
  if (ip.collect_dist && ip.distributions == nullptr)
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " ERROR: you need to provide a HistArray to compute distributions!" << std::endl;
      return 0;
    }
  // now everything is fine
  return 1;
}

Integrator::Integrator(std::ostream& ost):
  d_Ost(ost),
  d_GSLIntegrand{nullptr,0,nullptr},
  d_GSLRng(nullptr),
  d_GSLState(nullptr),
  d_externalGSLState(false)
{
}

Integrator::~Integrator()
{
  Reset();
}

void Integrator::Integrate(Integral& integral, integrand_par& ip, vegas_par& vp)
{
  if (!Init(integral,ip,vp))
    {
      d_Ost << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      d_Ost << " aborting integration ... " << std::endl;
      return;
    }

  int ncalls = vp.calls/vp.iterations;

  switch (ip.ps->whattype())
    {
    case 21:
      ncalls /= 10;
      if (vp.verbose>0) d_Ost << std::endl << "  PS 2->1 , #call/iter. = " << ncalls << std::endl;
      break;
    case 22:
      if (vp.verbose>0) d_Ost << std::endl << "  PS 2->2 , #call/iter. = " << ncalls << std::endl;
      break;
    case 23:
      ncalls *= 10;
      if (vp.verbose>0) d_Ost << std::endl << "  PS 2->3 , #call/iter. = " << ncalls << std::endl;
      break;
    default:
      if (vp.verbose>0) d_Ost << std::endl << "  unrecognized PS , #call/iter. = " << ncalls << std::endl;
    }

  int     ndim         = integral.GetDim();
  double* lower_limits = integral.Lo();
  double* upper_limits = integral.Up();
  std::vector<HistArray*>* dist   = ip.distributions;
  
  d_GSLState -> iterations = vp.iterations;
  d_GSLState -> verbose    = vp.verbose - 2;

  if (vp.grid_fixed)
    { // grid fixed -> no warmup run
      d_GSLState->alpha = 0;
    }
  else
    {
      /////////////////////////////////////////////////////////////////////////////////////////
      // WARMUP
      /////////////////////////////////////////////////////////////////////////////////////////
      double res = 0.0;
      double err = 0.0;
      // do not collect distributions in warmup run !
      if (dist != nullptr)
	{
	  for (auto e: *dist)
	    {
	      e->Pause();
	    }
	}
      
      d_GSLState->stage = 0; // reset grid, discard previous results
      gsl_monte_vegas_integrate (&d_GSLIntegrand, lower_limits, upper_limits, ndim, std::max(100000,ncalls), d_GSLRng, d_GSLState, &res, &err);
      if (vp.verbose >0)
	{
	  int prec = (int)log10(res/err)+1;
	  d_Ost << std::endl << " ===== Warmup-run ===== " << std::endl << "  res = " << std::setprecision(prec) << res << std::endl << "  err = " << std::setprecision(1) << (err) << std::endl;
	}
      /////////////////////////////////////////////////////////////////////////////////////////

    }

  //exit(1);
  /////////////////////////////////////////////////////////////////////////////////////////
  // MAIN RUN
  /////////////////////////////////////////////////////////////////////////////////////////
  d_GSLState->stage = 1; // re-use grid from previous run, only discard results
  d_GSLState->chisq = 0;
  // d_GSLState->alpha = 1.5;
  vp.result = 0;
  vp.error  = 0;

  if (vp.verbose >0) d_Ost << std::endl << " =====  Main-run  ===== " << std::endl;
  int run_count = 0;
  unsigned int INT_TIME     = 0;
  unsigned int INT_TIME_SUM = 0;

  // do max runs when collecting differential distributions...
  while ( (run_count<vp.max_runs) && (fabs(d_GSLState->chisq-1.0) > vp.chisq_limit) )
    {
      // activate histograms if available
      if (dist != nullptr)
	{
	  for (auto e: *dist)
	    {
	      e->Resume();
	    }
	}
      //////////// GSL-INTEGRATION ROUTINE ////////////
      INT_TIME      = clock(); // processor time
      gsl_monte_vegas_integrate (&d_GSLIntegrand, lower_limits, upper_limits, ndim, ncalls, d_GSLRng, d_GSLState,&vp.result, &vp.error);
      INT_TIME      = clock() - INT_TIME;
      INT_TIME_SUM += INT_TIME;
      /////////////////////////////////////////////////

      if (vp.verbose >0)
	{
	  int prec = (int)log10(vp.result/vp.error)+2;
	  d_Ost << "\n  (run " << run_count << ") ";
	  d_Ost << " i = " << std::setw(10) << std::setprecision(prec) << vp.result;
	  d_Ost << " s = " << std::setw(10) << std::setprecision(1) << (vp.error);
	  d_Ost << " (chisq/dof =  " << std::setw(4) << std::setprecision(2) << d_GSLState->chisq  << ") ";
	  d_Ost << " (time " << std::setw(4) << std::setprecision(1) << (float)INT_TIME/CLOCKS_PER_SEC << ") \n";
	}
      run_count++;
    }

  if (run_count == vp.max_runs)
    {
      if (vp.verbose >0) d_Ost << "\n  (max run limit reached, check results) \n";
    }

  // save run count for normalization of distributions
  vp.num_runs = run_count;
  vp.chisq    = d_GSLState->chisq;
  // normalize histograms to number of runs and iterations per run
  if (dist != nullptr)
    {
      for (auto e: *dist)
	{
	  // normalize distributions to #iterations, #runs (warm-up does not count)
	  e->Scale(vp.iterations*vp.num_runs);
	}
    }

  if (vp.verbose>-1)
    {
      int prec = (int)log10(vp.result/vp.error)+3;
      // print final result
      d_Ost << "\n ------------------------------------------------------------------";
      d_Ost << "\n | INTEGRAL = " << std::setw(11) << std::setprecision(prec) << vp.result << "  SIGMA = " << std::setw(11) << std::setprecision(2) << vp.error << " (chisq/dof = " << std::setw(4) << std::setprecision(2) << d_GSLState->chisq << ") |"; 
      d_Ost << "\n |  (total integration time " << std::setw(4) << std::setprecision(2) <<  (float)INT_TIME_SUM/CLOCKS_PER_SEC << " , mean " << std::setw(4) << std::setprecision(2) << (float)INT_TIME_SUM/CLOCKS_PER_SEC/run_count << ")                     |";
      d_Ost << "\n ------------------------------------------------------------------\n\n";
    }
  if (vp.verbose >0)
    {
      d_Ost << " ======================================================================\n\n";
    }
  // copy run parameters
  integral.LastRun = vp;
/////////////////////////////////////////////////////////////////////////////////////////
}
  

// set everything back to default values
void Integrator::Reset()
{
  d_GSLIntegrand = {nullptr,0,nullptr};
  if (d_GSLRng != nullptr)
    {
      //d_Ost << "\n  (deleting GSL rng " << d_GSLRng << ") \n";
      gsl_rng_free(d_GSLRng);
      d_GSLRng = nullptr;
    }
  if (!d_externalGSLState && d_GSLState != nullptr)
    {
      //d_Ost << "\n  (deleting GSL state " << d_GSLState << ") \n";
      gsl_monte_vegas_free(d_GSLState);
      d_GSLState = nullptr;
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
