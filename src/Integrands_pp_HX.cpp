
#include "../inc/Integrands_pp_HX.h"




///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

double Integrand_2_1_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  using namespace RunParameters;

  if (unlikely(dim!=2))
    {
      ERROR("wrong integral dimension!");
    }

  integrand_par* ip    = static_cast<integrand_par*> (arg);
  PS_2_1*        ps    = dynamic_cast<PS_2_1*> (ip->ps);
  unsigned long& flags = ip->eval_flags;  

  // hadronic c.m.e.
  double const& s_hadr = ip->s_hadr;
  // partonic c.m.e., fixed in this case
  double const& s_part = ps->get_msq(0);
  // double       rs_part = sqrt(s_part);
  // PDF integration
  double const& x0   = x[0];
  // the 2->1 phase space is proportional to delta(x0*x1*s_hadr-mH2),
  // x1-integration is eliminated
  double const& x1   = s_part/(s_hadr*x0);
  
  //ps->set();
  FV const& p1 = ps->p1();
  FV const& p2 = ps->p2();
  
  double res_b     = 0.0;
  double res_b_eff = 0.0;
  double res_v     = 0.0;
  double res_d     = 0.0;
  
  if (x1<1.0)
    {
      // if (x0<s_part/s_hadr || x1<s_part/s_hadr)
      // 	{
      // 	  PRINT(x0);
      // 	  PRINT(x1);
      // 	  PRINT(s_part/s_hadr);
      // 	}
      // extra factor 1/(x0*s_hadr) from phase space:
      // delta(s_part-mH^2) = 1/(x0*s_had)*delta(x1-mH^2/(x0*s_had))
      // using s_part = x0*x1*s_had
      double jac_flux_1 = TwoPi/(2.0*s_part*x0*s_hadr);

      // the PDF want dimensionful quantities in units of GeV (here, s_hadr is in units of mt!)
      double f1 = ip->pdf->xfxQ2(21, x0, MUF2*mScale2) / x0;
      double f2 = ip->pdf->xfxQ2(21, x1, MUF2*mScale2) / x1;
      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      double cf = PREF_GG*f1*f2*CONV_mt2i_pbarn;
      
      // Born
      if (flags & BOOST_BINARY(0 000 1))
	{
	  // PRINT(HiggsBosons::F_ggH1_s);
	  // PRINT(HiggsBosons::F_ggH1_s);
	  // PRINT(AmpPrefactors::PREF_B_PHIxQCD);
	  // full LO matrix elements
	  res_b     += Eval_B(p1,p2)*jac_flux_1*cf;
	  // PRINT(res_b);
	  // exit(1);
	  // effective LO matrix elements
	  res_b_eff += Eval_B_eff(p1,p2)*jac_flux_1*cf;
	}

      // virtual corrections
      if (flags & BOOST_BINARY(0 001 0))
      	{
      	  //finite part of the virtual corrections
      	  res_v += Eval_V(p1,p2)*jac_flux_1*cf;

      	  //finite part of the integrated dipoles
	  // (delta-terms minus +distribution endpoint terms)
      	  res_d += Eval_ID(p1,p2,x[1])*jac_flux_1*cf;
      	}


      // modified phase space for the integrated dipoles
      // is proportional to delta(x*x0*x1*s_hadr-mH2)
      // in this case we have 2.0*sp(p1_x,p2_x) = x*s_part = mH^2
      double const& x1_x = s_part/(s_hadr*x0*x[1]);
      //check if the modified ocnfiguration for integrated dipoles is physical
      if ( flags & BOOST_BINARY(0 001 0) && x1_x<1.0)
	{
	  // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
	  double jac_flux_x = TwoPi/(2.0*s_part*x[1]*x0*s_hadr);
	  double f1_x = ip->pdf->xfxQ2(21, x0  , MUF2*mScale2) / x0;
	  double f2_x = ip->pdf->xfxQ2(21, x1_x, MUF2*mScale2) / x1_x;

	  // include spin/color average of initial gluons and PFDs, convert to units of picobarn
	  double cf_x = PREF_GG*f1_x*f2_x*CONV_mt2i_pbarn;
	  res_d += Eval_ID_X(p1,p2,x[1])*jac_flux_x*cf_x;
	}
    }

  return (res_b + (ip->K)*(res_v + res_d));
}


///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

double Integrand_2_2_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  using namespace RunParameters;

  if (unlikely(dim!=3))
    {
      ERROR("wrong integral dimension!");
    }


  integrand_par* ip    = static_cast<integrand_par*> (arg);
  PS_2_2*        ps    = dynamic_cast<PS_2_2*> (ip->ps);
  unsigned long& flags = ip->eval_flags;


  // hadronic c.m.e.
  double s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];
  
  FV const& p1 = ps->p1();
  FV const& p2 = ps->p2();      
  FV const& p3 = ps->k2();// this should be the final state gluon -> m=0
      
  if (ps->set(sqrt(s_part),x[2],0.0))
    {
      // cuts on collinear and soft phase space regions (x[2] is the scattering angle)
      if (fabs(x[2]) > Cuts::COLL_CUT) return 0.0;
      if (p3[0]      < Cuts::SOFT_CUT) return 0.0;
      
      double res_r = 0.0;
      double res_d = 0.0;
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);
      
      // the PDF wants dimensionful quantities in units of GeV (here, s_hadr is in units of mt!)
      double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
      
      // real corrections gg - unitegrated dipoles
      if (flags & BOOST_BINARY(0 010 0))
	{
	  res_r += PREF_GG*f1_g*f2_g*jac_flux*(Eval_R_gg(p1,p2,p3)-Eval_UID(p1,p2,p3));
#ifdef DEBUG
	  CHECKNA(res_r);
#endif
	}
      // qq and qg contributions
      if (flags & BOOST_BINARY(1 100 0))
	{
	  static auto quark_pids     = { 5,  4,  3,  2,  1};
	  double qq_pdf = 0.0;
	  double qg_pdf = 0.0;
	  for (auto i: quark_pids)
	    {
	      double f1_q  = ip->pdf->xfxQ2(+i, x[0], MUF2*mScale2) / x[0];
	      double f2_qb = ip->pdf->xfxQ2(-i, x[1], MUF2*mScale2) / x[1];
	      // double f1_qb = ip->pdf->xfxQ2(-i, x[0], MUF2*mScale2) / x[0];
	      // double f2_q  = ip->pdf->xfxQ2(+i, x[1], MUF2*mScale2) / x[1];
	      qq_pdf += f1_q*f2_qb;//+f1_qb*f2_q;
	      qg_pdf += f1_q*f2_g+f1_g*f2_qb;
	    }

	  double z = ps->get_msq(0) / (2.0*sp(p1,p2));
	  if (flags & BOOST_BINARY(0 100 0))
	    {
	      //res_r += qq_pdf*PREF_QQ*Eval_R_qq(p1,p2,p3)*jac_flux;
	      res_r += qq_pdf*Eval_sigma_R_qq_fin(z);
#ifdef DEBUG
	      CHECKNA(res_r);
#endif
	    }
	  if (flags & BOOST_BINARY(1 000 0))
	    {
	      res_r += qg_pdf*Eval_sigma_R_qg_fin(z);
#ifdef DEBUG
	      CHECKNA(res_r);
#endif
	    }
	}
      // this will only be rescaled if the LO matrix elements were evaluated before !
      // otherwise K = 1
      return  (ip->K)*(res_r-res_d)*CONV_mt2i_pbarn;
    }
  else
    {
      return 0.0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////