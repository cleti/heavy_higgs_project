
#include "../inc/Integrands_pp_HX.h"





double Integrand_2_1_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  using namespace RunParameters;
  using namespace Bosons;
  
  integration_par* ip = static_cast<integration_par*> (arg);
  PS_2_1* ps          = dynamic_cast<PS_2_1*> (ip->ps);
  
  // hadronic c.m.e.
  double const& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double const& s_part = ps->get_msq();
  double       rs_part = sqrt(s_part);
  // PDF integration
  double const& x0   = x[0];
  // the 2->1 phase space is proportional to delta(x0*x1*s_hadr-mH2), x1-integration is eliminated
  double const& x1   = s_part/(s_hadr*x0);



  ps->set();
  FV const& p1 = ps->p1();
  FV const& p2 = ps->p2();
      
  double res_b = 0.0;
  double res_b_eff = 0.0;
  double res_v = 0.0;
  double res_d = 0.0;
	
  if (x1>0.0 && x1<1.0)
    {

      // factor 1/(x0*s_hadr) from phase space: delta(s_part-mH^2) = 1/(x0*s_had)*delta(x1-mH^2/(x0*s_had))
      double jac_flux_1 = 1.0/(2.0*s_part*x0*s_hadr);

      // the PDF want dimensionful quantities in units of GeV (here, s_hadr is in units of mt!)
      double f1 = ip->pdf->xfxQ2(21, x0, MUF2*mScale2) / x0;
      double f2 = ip->pdf->xfxQ2(21, x1, MUF2*mScale2) / x1;

      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      double cf = PREF_GG*f1*f2*CONV_mt2i_pbarn;

      if (CMP_LO)
	{
	  // full LO matrix elements
	  res_b     += Eval_B(p1,p2)*jac_flux_1*cf;
	  // effective LO matrix elements
	  res_b_eff += Eval_B_eff(p1,p2)*jac_flux_1*cf;
	}
      
      if (CMP_NLO_GG)
      	{
      	  //finite part of the virtual corrections
      	  res_v += Eval_V(p1,p2)*jac_flux_1*cf;

      	  //finite part of the integrated dipoles (delta-terms minus +distribution endpoint terms)
      	  res_d += Eval_ID(p1,p2,x[1])*jac_flux_1*cf;
      	}
    }

  // modified phase space for the integrated dipoles is proportional to delta(x*x0*x1*s_hadr-mH2)
  // in this case we have 2.0*sp(p1_x,p2_x) = x*s_part = mH^2
  double const& x1_x = M2_1/(s_hadr*x0*x[1]);
  //check if the modified ocnfiguration for integrated dipoles is physical
  if (CMP_NLO_GG && x1_x>0.0 && x1_x<1.0)
    {
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_x = 1.0/(2.0*s_part*x[1]*x0*s_hadr);
      double f1_x = ip->pdf->xfxQ2(21, x0  , MUF2*mScale2) / x0;
      double f2_x = ip->pdf->xfxQ2(21, x1_x, MUF2*mScale2) / x1_x;

      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      double cf_x = PREF_GG*f1_x*f2_x*CONV_mt2i_pbarn;
      res_d += Eval_ID_X(p1,p2,x[1])*jac_flux_x*cf_x;
    }

  return (res_b) + (res_v + res_d);
}
