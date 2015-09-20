

#include "../inc/Integrands_pp_ttX.h"




// these numbers indentify the quark PDF sets in LHAPDF
// 5: b-quark ... 1: d-quark, negative numbers -> antiquarks
const int quark_pids[] = {  5,  4,  3,  2,  1};

  


#ifdef WITH_T_SPIN
static int Eval_tt_decay(
			 PS_2_2* ps,
			 const HiggsModel& hm,
			 const double* x,
			 double& tt_decay)
{
  using namespace Constants;
  tt_decay = 1.0;
  /////////////////////////////////////////////////////////////////////////////////////
  // polarized amplitudes: evaluate decay phase spaces and t/tbar decay matrix elements
  /////////////////////////////////////////////////////////////////////////////////////
  double const& mt2    = hm.mt2();

  static PS_2_3* ps_t1 = static_cast<PS_2_3*>(ps->get_child(0));
  static PS_2_3* ps_t2 = static_cast<PS_2_3*>(ps->get_child(1));

#ifdef DEBUG
  if ( ps_t1 == nullptr ) {ERROR("'ps_t1': nullptr received");}
  if ( ps_t2 == nullptr ) {ERROR("'ps_t2': nullptr received");}
#endif
  
  // the function Eval_t_blnu() will store in S1_r, S2_r
  // the lepton momenta defined in the t/tbar restframes 
  FV& S1_r = ps->s1_r();
  FV& S2_r = ps->s2_r();

  if ( ps_t1->set(mt2,x[3] ,x[4] ,x[5] ,x[6] ,x[7] ) &&
       ps_t2->set(mt2,x[8] ,x[9] ,x[10],x[11],x[12])    )
    {
      ///////////////// decay of top quark ////////////////////////////////////////////
      tt_decay *= (ps_t1->get_wgt() * Eval_t_blnu(*ps_t1,S1_r,hm));
      ///////////////// decay of anti-top quark ///////////////////////////////////////
      tt_decay *= (ps_t2->get_wgt() * Eval_t_blnu(*ps_t2,S2_r,hm));
      /////////////////////////////////////////////////////////////////////////////////
      S1_r *=  kappa_p;
      S2_r *=  kappa_m;// antitop spin vector receives relative sign through kappa_m!
      /////////////////////////////////////////////////////////////////////////////////
      // in the production matrix element we have to subtitute S1 -> kappa_p * lepton-mom.
      // and  S1 -> kappa_m * antilepton-mom., but these factors have to be removed again
      // when observables involving the lepton/antilepton momenta are considered !!!
      /////////////////////////////////////////////////////////////////////////////////
      // normalization factor from t/tbar propagators pi/(mt*GammaT) * 1/(2*pi)
      // factor 1/GammaT already included in the output of Eval_t_blnu() !
      tt_decay /= (4.0*mt2); 
    }
  else
    {
      return 0;
    }

  ////////////////////////////////////////////////////////////////////
  // boost top spin and decay vectors to parton z.m.f. ///////////////
  ////////////////////////////////////////////////////////////////////
  // top anti-top 4-vectors given in the tt z.m.f.
  FV const& k1 = ps->k1();
  FV const& k2 = ps->k2();
  FV& S1 = ps->s1();
  FV& S2 = ps->s2();
  // FV K1 = {1,0,0,0};
  // FV K2 = {1,0,0,0};
  // now we have to boost the spin vectors to the parton z.m.f. = ttbar z.m.f. 
  S1 = S1_r;
  S2 = S2_r;
  // boost from k1 restframe to tt z.m.f.
  static LT boost_k1_RF;
  boost_k1_RF.set_boost(k1,1);
  // boost from k2 restframe to tt z.m.f.
  static LT boost_k2_RF;
  boost_k2_RF.set_boost(k2,1);
  
  boost_k1_RF.apply(S1);
  boost_k2_RF.apply(S2);

  
  // PRINT(sp(k1,S1));
  // PRINT(sp(k2,S2));
  // PRINT_4VEC(S1);
  // PRINT_4VEC(S2);
  
  // FV S1t;
  // FV S2t;
  // set_spins_in_tt_zmf(k1,k2,S1t,S2t,S1_r,S2_r);
  // PRINT_4VEC(S1_r);
  // PRINT_4VEC(S2_r);
  // PRINT_4VEC(S1t);
  // PRINT_4VEC(S2t);
  // PRINT_4VEC(S1);
  // PRINT_4VEC(S2);
  // exit(1);
  ////////////////////////////////////////////////////////////////////
  tt_decay *= 4.0;//*0.7983;
  return tt_decay;
}
#endif





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_poly2(double* x, size_t dim, void* arg)
{
  integrand_par* ip    = static_cast<integrand_par*> (arg);
  DistVec* dist = (ip->distributions);
  double vwgt  = ip->cmp_v_weight();
  double res = 1.0+x[0]+x[1]*x[1];
  (*dist)[0]->GetHistograms()->FillOne(H_LO_QCD,x[1],res*vwgt);
  return res;
}

// to reproduce Bernies s_part plots
double Integrand_2_2(double* x, size_t dim, void* arg)
{
  using namespace Constants;

  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_2*         ps    = dynamic_cast<PS_2_2*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   mScale = hm->Scale();  
  double const&   mScale2= hm->Scale2();

  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;
  
  double rs_part = 0.0;
  double y = 0.0;
  if (dim==1)
    {
      rs_part = sqrt(ip->s_hadr);
      y = x[0];
    }
  else if (dim==2)
    {
      rs_part = x[0];
      y = x[1];
    }
  else
    {
      ERROR("wrong integral dimension!");
    }

  if (ps->set(rs_part,y,0.0))
    {
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*rs_part*rs_part);

      // born matrix elements (GG)
      double res_b = Eval_B(*ps,*hm,flags,0)*jac_flux*cf_gg;
      
      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  auto dist = (ip->distributions);
	  FV const& p1 = ps->p1();
	  FV const& p2 = ps->p2();	  
	  double vwgt  = ip->cmp_v_weight();
	  if (flags & F_EVAL_B_QCDxQCD)
	    { // only the s_part distribution
	      (*dist)[0]->GetHistograms()->FillOne(H_LO_QCD,obs_M12(p1,p2)*mScale,res_b*vwgt);
	    }
	  else
	    { // only the s_part distribution 
	      (*dist)[0]->GetHistograms()->FillOne(H_LO_PHI,obs_M12(p1,p2)*mScale,res_b*vwgt);
	    }
	}
      ////////////////////////////////////////////////////////////////////////////
      return  (res_b);
    }
  else
    {
      return 0.0;
    }
  
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_2_2_pdf_BV(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef DEBUG
#ifdef WITH_T_SPIN
    if ( dim!=13 )
#else
    if ( dim!=3 )
#endif
      {
	ERROR("wrong integral dimension!");
      }
#endif
  

  // some aliases
  integrand_par*    ip      = static_cast<integrand_par*> (arg);
  double const      &s_hadr = ip->s_hadr;
  PS_2_2*const      &ps     = ip->ps_2;
#ifdef DEBUG
  if ( ps == nullptr )
    {
      ERROR("ps_2 = nullptr in integrand parameters!");
    }
#endif  
  ulong const       &flags  = ip->eval_flags;
  HiggsModel*const  &hm     = ip->higgs_model;
  CutVec*const      &cuts   = ip->cuts;
  
  double const      &MUF2   = hm->MUF2();
  double const      &mScale = hm->Scale(); 
  double const      &mScale2= hm->Scale2();
    
  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];

  // set up 2->2 phase space, return 0.0 for unphysical settings (s < 4 mt^2)
  if (!(ps->set(sqrt(s_part),x[2],0.0))) return 0.0;
  	    
  // evaluate phase space cuts, ps is defined in the parton z.m.f. at this point
  if (!EvalCuts(cuts,ps,nullptr)) return 0.0;

  // need this only in case of spin dependent amplitudes
  double tt_decay = 1.0;
#ifdef WITH_T_SPIN
  if (!Eval_tt_decay(ps,*hm,x,tt_decay)) exit(1);
#endif
  
  // use eff. Higgs-gluon coupling in the ME if rescaling factor is used
  bool EFF = false; // hm->UseK();
  
  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;
  
  // spin/color average for qq initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_qq = PREF_QQ*CONV_GeV2i_pbarn/mScale2;

  // results
  double res_b = 0.0;
  double res_v = 0.0;
  // needed for distributions
  double   vwgt = 0.0;
  DistVec* dist = nullptr;
  if (ip->collect_dist)
    {
      dist = ip->distributions;
      vwgt = ip->cmp_v_weight();
      // note: these are the proton momenta in the parton z.m.f.!!!
      // needed to compute the boost to the lab frame
      ps->P1() = 1.0/x[0] * (ps->p1());
      ps->P2() = 1.0/x[1] * (ps->p2());
    }
  
  // gluon PDFs
  // scale MU must be given in [GeV] !!! 
  double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
  double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

  // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
  double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);

  // born matrix elements (gg initial state)
  if (EVAL_B(flags))
    {
      res_b += Eval_B(*ps,*hm,flags,EFF)*jac_flux*cf_gg*f1_g*f2_g*tt_decay;
    }
       
  // born matrix elements (qq initial state)
  if (flags & F_EVAL_B_QCDxQCD)
    {  
      double qq_pdf = 0.0;
      // get the quark/antiquark pdfs
      for (auto i: quark_pids)
	{
	  double f1_q  = ip->pdf->xfxQ2(+i, x[0], MUF2*mScale2) / x[0];
	  double f2_qb = ip->pdf->xfxQ2(-i, x[1], MUF2*mScale2) / x[1];
	  double f1_qb = ip->pdf->xfxQ2(-i, x[0], MUF2*mScale2) / x[0];
	  double f2_q  = ip->pdf->xfxQ2(+i, x[1], MUF2*mScale2) / x[1];
	  qq_pdf += (f1_q*f2_qb+f1_qb*f2_q);
	}
      // the LO QCD Born amplitudes are symmetric under y <-> -y
      // (also with spin dependence!)
      // --> PDFs can be summed and pulled out
      res_b += Eval_B_QQ(*ps,*hm)*jac_flux*cf_qq*qq_pdf*tt_decay;
    }

  if (ip->collect_dist)
    {
      // we put the result into the LO_QCD histogram if the according
      // flag is set, in the LO_PHI histogram otherwise
      H_Index ih = (flags & F_EVAL_B_QCDxQCD)?H_LO_QCD:H_LO_PHI;
      ps->FillDistributions(*dist,ih,res_b*vwgt,mScale);
    }
  
  // finite part of the virtual corrections (gg initial state)
  if (EVAL_V(flags))
    {
      res_v += Eval_V(*ps,*hm,flags)*jac_flux*cf_gg*f1_g*f2_g*tt_decay;
      
      if (ip->collect_dist)
	{
	  ps->FillDistributions(*dist,H_NLO_PHI_V,res_v*vwgt,mScale);
	}
    }

  return (res_b+res_v);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_2_2_pdf_ID(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef DEBUG
#ifdef WITH_T_SPIN
  if ( dim!=14 )
#else
    if ( dim!=4 )
#endif
      {
	ERROR("wrong integral dimension!");
      }
#endif

  // some aliases
  integrand_par*    ip      = static_cast<integrand_par*> (arg);
  double const      &s_hadr = ip->s_hadr;
  PS_2_2*const      &ps     = ip->ps_2;
#ifdef DEBUG
  if ( ps == nullptr )
    {
      ERROR("ps_2 = nullptr in integrand parameters!");
    }
#endif  
  ulong const       &flags  = ip->eval_flags;
  HiggsModel*const  &hm     = ip->higgs_model;
  CutVec*const      &cuts   = ip->cuts;
  
  double const      &MUF2   = hm->MUF2();
  double const      &mt2    = hm->mt2();
  double const      &mScale = hm->Scale(); 
  double const      &mScale2= hm->Scale2();


  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];

  // set up the unmodified 2->2 phase space, return 0.0 for unphysical settings (s < 4 mt^2)
  // since x < 1 => s' = x*s < s, thus the setting of the modified phase spaces
  // can only be physical if that of the unmodified already is
  if (! ps->set(sqrt(s_part),x[2],0.0)) return 0.0;
  
  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;

  // results
  double res_d = 0.0;
  double res_dx1 = 0.0;
  double res_dx2 = 0.0;
  // needed for distributions
  double   vwgt = 0.0;
  DistVec* dist = nullptr;
  if (ip->collect_dist)
    {
      dist = ip->distributions;
      vwgt = ip->cmp_v_weight();
      // note: these are the proton momenta in the parton z.m.f.!!!
      // needed to compute the boost to the lab frame
      ps->P1() = 1.0/x[0] * (ps->p1());
      ps->P2() = 1.0/x[1] * (ps->p2());
    }
  
  // gluon PDFs
  // scale MU must be given in [GeV] !!! 
  double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
  double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

  // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
  double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*s_part);
      
  // finite part of the integrated dipoles: delta terms and distribution end-point terms
  // evaluate phase space cuts, ps is defined in the parton z.m.f. at this point
  if (EvalCuts(cuts,ps,nullptr))
    {
      res_d += Eval_ID_GG(*ps,x[3],*hm,flags)*jac_flux_1*cf_gg*f1_g*f2_g;
      
      if (ip->collect_dist)
	{
	  ps->FillDistributions(*dist,H_NLO_PHI_ID,(res_d)*vwgt,mScale);
	}
    }



  /*
    the continuum part of the integrated dipoles is evaluated on
    a boosted phase space, the boost is along the +-z directions,
    x[3]=(1+beta)/(1-beta) is the boost parameter, the parton c.m.e.
    is then x*s
  */
      
  // boosted phase space (boost is done in Eval_ID_X)
  static PS_2_2 ps_x1(mt2,mt2,"x*p1 p2 -> k1' k2'  [static in Integrand_2_2_pdf_ID]");
  static PS_2_2 ps_x2(mt2,mt2,"p1 x*p2 -> k1' k2'  [static in Integrand_2_2_pdf_ID]");
  static LT dip_boost;

  // finite part of the integrated dipoles: x-dependent terms of P and K operators
  double sx = x[3]*s_part;
  // first set up boosted phase spaces in the x*p1-p2 / p1-x*p2 z.m.f.
  if ( ps_x1.set(sqrt(sx),x[2],0.0) && x[3]<Cuts::IDIP_X_CUT )
    {
      // modified version of flux and phase space density
      double jac_flux_x = TwoPi*ps_x1.get_wgt()/(2.0*sx);
      
      // the c.m.e. s'=x*s is equal for ps_x1 and ps_x2,
      // thus the 4 vectors in the x*p1-p2 / p1-x*p2 z.m.f. are also the same
      ps_x2 = ps_x1;

      // ps_x1.print();
      // ps_x2.print();
      
      // now we boost the 4-vectors to the p1-p2 z.m.f. (the frame in which ps is set up)
      // boost in +z direction is applied to ps_x1 4-vectors (x*p1-p2 z.m.f. -> p1-p2 z.m.f.)
      dip_boost.set_boost_z(x[3],0);
      dip_boost.apply(ps_x1.p1());
      dip_boost.apply(ps_x1.p2());
      dip_boost.apply(ps_x1.k1());
      dip_boost.apply(ps_x1.k2());
#ifdef WITH_T_SPIN
      dip_boost.apply(ps_x1.s1());
      dip_boost.apply(ps_x1.s2());
#endif

      // boost in -z direction is applied to ps_x2 4-vectors (p1-x*p2 z.m.f. -> p1-p2 z.m.f.)
      dip_boost.set_boost_z(x[3],1);
      dip_boost.apply(ps_x2.p1());
      dip_boost.apply(ps_x2.p2());
      dip_boost.apply(ps_x2.k1());
      dip_boost.apply(ps_x2.k2());
#ifdef WITH_T_SPIN
      dip_boost.apply(ps_x2.s1());
      dip_boost.apply(ps_x2.s2());
#endif

      // ps_x1.print();
      // ps_x2.print();

      // EXIT(1);
      
      // evaluate cuts on the modified phase spaces,
      // ps_x1 and ps_x2 are defined in the parton z.m.f. at this point
      bool cuts_x1 = EvalCuts(cuts,&ps_x1,nullptr);
      bool cuts_x2 = EvalCuts(cuts,&ps_x2,nullptr);

      // evaluate gg contribution on both phase spaces
      // this is always symmetric under p1 <-> p2
      if (flags & F_EVAL_D_GG_CONT)
	{
	  if (cuts_x1)
	    {
	      res_dx1 += 0.5*Eval_ID_X_GG(ps_x1,x[3],*hm,flags)*jac_flux_x*cf_gg*f1_g*f2_g;
	    }
	  if (cuts_x2)
	    {
	      res_dx2 += 0.5*Eval_ID_X_GG(ps_x2,x[3],*hm,flags)*jac_flux_x*cf_gg*f1_g*f2_g;
	    }
	}

      // evaluate qg/gq contribution on both phase spaces
      // not symmetric under p1 <-> p2 if spin dependence is included
      if (flags & F_EVAL_D_QG_CONT)
	{	  
	  double qg_pdf = 0.0;
	  double gq_pdf = 0.0;
	  for (auto i: quark_pids)
	    {
	      double f1_q  = ip->pdf->xfxQ2(+i, x[0], MUF2*mScale2) / x[0];
	      double f2_qb = ip->pdf->xfxQ2(-i, x[1], MUF2*mScale2) / x[1];
	      double f1_qb = ip->pdf->xfxQ2(-i, x[0], MUF2*mScale2) / x[0];
	      double f2_q  = ip->pdf->xfxQ2(+i, x[1], MUF2*mScale2) / x[1];
	      qg_pdf += f2_g*(f1_q+f1_qb);
	      gq_pdf += f1_g*(f2_q+f2_qb);
	    }


	  if (cuts_x1) // this assumes that the cuts are invariant under p1 <-> p2
	    {
#ifdef WITH_T_SPIN
	      // the spin-dependent matrix elements are not symmetric under p1 <-> p2
	      res_dx1 += 0.5*Eval_ID_X_QG(ps_x1,x[3],*hm,flags)*jac_flux_x*cf_gg*qg_pdf;
	      ps_x1.swap_initial_state(); // p1 <-> p2
	      res_dx1 += 0.5*Eval_ID_X_QG(ps_x1,x[3],*hm,flags)*jac_flux_x*cf_gg*gq_pdf;
	      ps_x1.swap_initial_state(); // p2 <-> p1
#else
	      res_dx1 += 0.5*Eval_ID_X_QG(ps_x1,x[3],*hm,flags)*jac_flux_x*cf_gg*(qg_pdf+gq_pdf);
#endif
	    }

	  if (cuts_x2) // this assumes that the cuts are invariant under p1 <-> p2
	    {
#ifdef WITH_T_SPIN
	      // the spin-dependent matrix elements are not symmetric under p1 <-> p2
	      res_dx2 += 0.5*Eval_ID_X_QG(ps_x2,x[3],*hm,flags)*jac_flux_x*cf_gg*qg_pdf;
	      ps_x2.swap_initial_state(); // p1 <-> p2
	      res_dx2 += 0.5*Eval_ID_X_QG(ps_x2,x[3],*hm,flags)*jac_flux_x*cf_gg*gq_pdf;
	      ps_x2.swap_initial_state(); // p2 <-> p1
#else
	      res_dx2 += 0.5*Eval_ID_X_QG(ps_x2,x[3],*hm,flags)*jac_flux_x*cf_gg*(qg_pdf+gq_pdf);
#endif
	    }
	}
	  
      // DISTRIBUTIONS: x neq 1 contribution //////////////////////////////////////
      if (ip->collect_dist)
	{
	  // at this point, 4-vectors in ps_x1 / ps_x2 should be
	  // in the same reference frame as ps
	  // -> can copy proton momenta from ps
	  ps_x1.P1() = (ps->P1());
	  ps_x1.P2() = (ps->P2());
	  ps_x2.P1() = (ps->P1());
	  ps_x2.P2() = (ps->P2());
	  // ps_x1.set_x(x[3]);
	  // ps_x2 = ps_x1;
	      
	  // dip_boost.set_boost_z(x[3],0);
	  // dip_boost.apply(ps_x1.k1());
	  // dip_boost.apply(ps_x1.k2());
	      
	  // dip_boost.set_boost_z(x[3],1);
	  // dip_boost.apply(ps_x2.k1());
	  // dip_boost.apply(ps_x2.k2());
	  if (cuts_x1)
	    {
	      ps_x1.FillDistributions(*dist,H_NLO_PHI_ID,(res_dx1)*vwgt,mScale);
	    }

	  if (cuts_x2)
	    {
	      ps_x2.FillDistributions(*dist,H_NLO_PHI_ID,(res_dx2)*vwgt,mScale);
	    }
	}
      ////////////////////////////////////////////////////////////////////////////
    }
    
  return  (res_d+res_dx1+res_dx2);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_2_3_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef DEBUG
#ifdef WITH_T_SPIN
  if ( dim !=16 )
#else
    if ( dim != 6 )
#endif
      {
	ERROR("wrong integral dimension!");
      }
#endif
  
  integrand_par* ip     = static_cast<integrand_par*> (arg);
  double const&  s_hadr = ip->s_hadr;
  PS_2_3*&       ps     = ip->ps_3;
#ifdef DEBUG
  if ( ps == nullptr )
    {
      ERROR("ps_3 = nullptr in integrand parameters!");
    }
#endif  
  const ulong&   flags  = ip->eval_flags;
  HiggsModel*&   hm     = ip->higgs_model;
  CutVec*&       cuts   = ip->cuts;
  
  double const&  MUF2   = hm->MUF2();
  double const&  mScale = hm->Scale(); 
  double const&  mScale2= hm->Scale2();
  
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];
  
  // set up the 2->3 phase space, return 0.0 for unphysical settings
  if (! ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5])) return 0.0;

  // technical cut on soft/collinear phase space region
  if (fabs(x[2])  > Cuts::COLL_CUT) return 0.0;
  if (ps->p3()[0] < Cuts::SOFT_CUT) return 0.0;
  
  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;

  // results
  double res_r_gg = 0.0;
  double res_d_gg = 0.0;
  // needed for distributions
  double   vwgt = 0.0;
  DistVec* dist = nullptr;
  if (ip->collect_dist)
    {
      dist = ip->distributions;
      vwgt = ip->cmp_v_weight();
      // note: these are the proton momenta in the parton z.m.f.!!!
      // needed to compute the boost to the lab frame
      ps->P1() = 1.0/x[0] * (ps->p1());
      ps->P2() = 1.0/x[1] * (ps->p2());
    }

  // the PDF wants dimensionful quantities in units of GeV (here, MUF is in units of mt!) 
  double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
  double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

  // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
  double jac_flux_gg = TwoPi*ps->get_wgt()/(2.0*s_part)*cf_gg*f1_g*f2_g;
      
  // evaluate the phase space cuts and add the real contribution if passed
  // ps is defined in the parton z.m.f. at this point
  if (EvalCuts(cuts,ps,nullptr))
    {
      res_r_gg += Eval_R_GG(*ps,*hm,flags)*jac_flux_gg;
      
      if (ip->collect_dist)
	{
	  ps->FillDistributions(*dist,H_NLO_PHI_R,res_r_gg*vwgt,mScale);
	}
    }

  
  // the contribution of the unintegrated dipoles (cuts checked for each dipole config.)
  res_d_gg += Eval_UID_GG (*ps,*hm,flags,-vwgt*jac_flux_gg,dist,cuts)*jac_flux_gg;
  
  if ( ps->p3()[0] < 1e-3 || (1.0-fabs(x[2])) < 1e-3 )
    {
      // if (fabs((res_d_gg-res_r_gg)/(res_d_gg+res_r_gg)) > 1e-5)
      // 	{
      // 	  PRINT((1.0-fabs(x[2])));
      // 	  PRINT(ps->p3()[0]);
      // 	  PRINT(res_r_gg);
      // 	  PRINT(res_d_gg);
      // 	}
      g_hist_r_wgts.Fill(std::log(std::fabs(vwgt*res_r_gg)));
      g_hist_uid_wgts.Fill(std::log(std::fabs(vwgt*res_d_gg)));
    }
  
  return (res_r_gg-res_d_gg);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_2_3_qg_qq_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef DEBUG
#ifdef WITH_T_SPIN
  if ( dim != 16 )
#else
    if ( dim != 6 )
#endif
      {
	ERROR("wrong integral dimension!");
      }
#endif
  
  integrand_par* ip     = static_cast<integrand_par*> (arg);
  double const&  s_hadr = ip->s_hadr;
  PS_2_3*&       ps     = ip->ps_3;
#ifdef DEBUG
  if ( ps == nullptr )
    {
      ERROR("ps_3 = nullptr in integrand parameters!");
    }
#endif
  ulong const&   flags  = ip->eval_flags;
  HiggsModel*&   hm     = ip->higgs_model;
  CutVec*const&  cuts   = ip->cuts;

  double const&  MUF2   = hm->MUF2();
  double const&  mScale = hm->Scale();
  double const&  mScale2= hm->Scale2();

  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];

  // set up the 2->3 phase space, return 0.0 for unphysical settings
  if (!ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5])) return 0.0;

  // technical cuts on soft/collinear phase space
  if (fabs(x[2])  > Cuts::COLL_CUT) return 0.0;
  if (ps->p3()[0] < Cuts::SOFT_CUT) return 0.0;

  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;

  // spin/color average for qg/gq initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_qg = PREF_QG*CONV_GeV2i_pbarn/mScale2;

  // spin/color average for qq initial state, conversion GeV^-2 -> picobarn
  // this assumes that mScale does not change during integration!
  static const double cf_qq = PREF_QQ*CONV_GeV2i_pbarn/mScale2;

  // results
  double res_r_qq = 0.0;
  double res_r_qg = 0.0;
  double res_d_qg = 0.0;
  // needed for distributions
  double   vwgt = 0.0;
  DistVec* dist = nullptr;
  if (ip->collect_dist)
    {
      dist = ip->distributions;
      vwgt = ip->cmp_v_weight();
      // note: these are the proton momenta in the parton z.m.f.!!!
      // needed to compute the boost to the lab frame
      ps->P1() = 1.0/x[0] * (ps->p1());
      ps->P2() = 1.0/x[1] * (ps->p2());
    }
      
  // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
  double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);

  /////////////////////////// PDF //////////////////////////////////
  double qq_pdf = 0.0;
  double qg_pdf = 0.0;
  double gq_pdf = 0.0;
  double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
  double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
  // get the quark/antiquark pdfs
  for (auto i: quark_pids)
    {
      double f1_q  = ip->pdf->xfxQ2(+i, x[0], MUF2*mScale2) / x[0];
      double f2_qb = ip->pdf->xfxQ2(-i, x[1], MUF2*mScale2) / x[1];
      double f1_qb = ip->pdf->xfxQ2(-i, x[0], MUF2*mScale2) / x[0];
      double f2_q  = ip->pdf->xfxQ2(+i, x[1], MUF2*mScale2) / x[1];
      qq_pdf += f1_q*f2_qb+f1_qb*f2_q;
      qg_pdf += f2_g*(f1_q+f1_qb);
      gq_pdf += f1_g*(f2_q+f2_qb);
    }
  //////////////////////////////////////////////////////////////////

  // evaluate phase space cuts, ps is defined in the parton z.m.f. at this point
  if ( EvalCuts(cuts,ps,nullptr) )
    {
      res_r_qq += Eval_R_QQ(*ps,*hm,flags)*cf_qq*jac_flux*qq_pdf;
      res_r_qg += Eval_R_QG(*ps,*hm,flags)*cf_qg*jac_flux*qg_pdf;

      // crossed configuration
      ps->swap_initial_state(); // p1 <-> p2
      res_r_qg += Eval_R_QG(*ps,*hm,flags)*cf_qg*jac_flux*gq_pdf;
      ps->swap_initial_state(); // p2 <-> p1
	  
      if (ip->collect_dist)
	{
	  ps->FillDistributions(*dist,H_NLO_PHI_R,(res_r_qq+res_r_qg)*vwgt,mScale);
	}
    }
      
  // the contribution of the unintegrated qg dipoles (cuts checked for each dipole config.)
  res_d_qg += Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_gg*jac_flux*qg_pdf,dist,cuts)*cf_gg*jac_flux*qg_pdf;
  ps->swap_initial_state(); // p1 <-> p2
  res_d_qg += Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_gg*jac_flux*gq_pdf,dist,cuts)*cf_gg*jac_flux*gq_pdf;
  ps->swap_initial_state(); // p2 <-> p1
      
  return res_r_qq + (res_r_qg-res_d_qg);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
