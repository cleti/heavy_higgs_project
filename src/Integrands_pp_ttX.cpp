

#include "../inc/Integrands_pp_ttX.h"



static int Eval_tt_decay(
			 PS_2_2* ps,
			 const HiggsModel& hm,
			 const double* x,
			 double& tt_decay)
{
  using namespace Constants;
  tt_decay = 1.0;
#ifdef WITH_T_SPIN
  /////////////////////////////////////////////////////////////////////////////////////
  // polarized amplitudes: evaluate decay phase spaces and t/tbar decay matrix elements
  /////////////////////////////////////////////////////////////////////////////////////
  double const& mt2    = hm.mt2();

  static PS_2_3* ps_t1 = static_cast<PS_2_3*>(ps->get_child(0));
  static PS_2_3* ps_t2 = static_cast<PS_2_3*>(ps->get_child(1));
	  
  if (unlikely(ps_t1 == nullptr)) {ERROR("'ps_t1': nullptr received");}
  if (unlikely(ps_t2 == nullptr)) {ERROR("'ps_t2': nullptr received");}

  // the function Eval_t_blnu() will store in S1_r, S2_r
  // the lepton momenta defined in the t/tbar restframes 
  FV& S1_r = ps->s1_r();
  FV& S2_r = ps->s2_r();


  if ( ps_t1->set(mt2,x[3] ,x[4] ,x[5] ,x[6] ,x[7] ) &&
       ps_t2->set(mt2,x[8] ,x[9] ,x[10],x[11],x[12])    )
    {
      // factor 1/GammaT already included in the output of Eval_t_blnu() !
      ///////////////// decay of top quark ////////////////////////////////////////////
      tt_decay *= (ps_t1->get_wgt() * Eval_t_blnu(*ps_t1,S1_r,hm));
      // spin vector in top restframe
      /////////////////////////////////////////////////////////////////////////////////
      //PRINT(ps_t1->get_wgt());PRINT(Eval_t_blnu(*ps_t1,S1_r,hm));
      ///////////////// decay of anti-top quark ///////////////////////////////////////
      tt_decay *= (ps_t2->get_wgt() * Eval_t_blnu(*ps_t2,S2_r,hm));
      // PRINT_4VEC(S2_r);
      S2_r *= (-1);// antitop spin vector receives relative sign
      // // PRINT_4VEC(S2_r);
      // // exit(1);
      /////////////////////////////////////////////////////////////////////////////////
      //PRINT(ps_t2->get_wgt());PRINT(Eval_t_blnu(*ps_t2,S1_r,hm));
      // normalization factor from t/tbar propagators pi/(mt*GammaT) * 1/(2*pi)
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
  // S2 *= -1;
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
#else
  /////////////////////////////////////////////////////////////////////////////////////
  // unpolarized amplitudes: mulitply tt -> l + jets branching fraction, approx. 24/81
  /////////////////////////////////////////////////////////////////////////////////////
  // factor 4 accounts for summation over tt spin configurations
  //tt_decay *= 24./81.*4.0;
#endif
  return 1;
}





//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_poly2(double* x, size_t dim, void* arg)
{
  integrand_par* ip    = static_cast<integrand_par*> (arg);
  std::vector<HistArray*>* dist = (ip->distributions);
  double vwgt  = ip->cmp_v_weight();
  double res = 1.0+x[0]+x[1]*x[1];
  (*dist)[0]->FillOne(H_LO_QCD,x[1],res*vwgt);
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
	      (*dist)[0]->FillOne(H_LO_QCD,obs_M12(p1,p2)*mScale,res_b*vwgt);
	    }
	  else
	    { // only the s_part distribution 
	      (*dist)[0]->FillOne(H_LO_PHI,obs_M12(p1,p2)*mScale,res_b*vwgt);
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
double Integrand_2_2_pdf_BV(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef WITH_T_SPIN
  if (unlikely(dim!=13))
#else
  if (unlikely(dim!=3))
#endif
  {
    ERROR("wrong integral dimension!");
  }

  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_2*         ps    = dynamic_cast<PS_2_2*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   MUF2   = hm->MUF2();
  double const&   mScale = hm->Scale(); 
  double const&   mScale2= hm->Scale2();
  // use eff. Higgs-gluon coupling in the ME if rescaling factor is used
  bool EFF = hm->UseK();
  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;
  // spin/color average for qq initial state, conversion GeV^-2 -> picobarn
  static const double cf_qq = PREF_QQ*CONV_GeV2i_pbarn/mScale2;
    
  // hadronic c.m.e.
  double& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];

  /////////////////////////////////////////////////////
  // Mtt cut //////////////////////////////////////////
  /////////////////////////////////////////////////////
#ifdef WITH_T_SPIN
  // double const& MH = hm->GetBoson(0)->M();
  // if (fabs(sqrt(s_part)-MH)>10.0/mScale) return 0.0;
#endif
  /////////////////////////////////////////////////////
  
  if (ps->set(sqrt(s_part),x[2],0.0))
    {
      double tt_decay = 1.0;
      if (!Eval_tt_decay(ps,*hm,x,tt_decay)) exit(1);
      
      double res_b = 0.0;
      double res_v = 0.0;

      // gluon PDFs
      // scale MU must be given in [GeV] !!! 
      double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);
      
      // born matrix elements (GG)
      if (EVAL_B(flags)) res_b = Eval_B(*ps,*hm,flags,EFF)*jac_flux*cf_gg*f1_g*f2_g*tt_decay;

      // born matrix elements (QQ)
      if (flags & F_EVAL_B_QCDxQCD)
        {
          static auto quark_pids     = {  5,  4,  3,  2,  1};	  
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
          // the QCD Born amplitudes are symmetric under y <> -y
          // --> PDFs can be pulled out
          res_b += Eval_B_QQ(*ps,*hm)*jac_flux*cf_qq*qq_pdf*tt_decay; 
        }
      // finite part of the virtual corrections (GG)
      if (EVAL_V(flags)) res_v = Eval_V(*ps,*hm,flags)*jac_flux*cf_gg*f1_g*f2_g*tt_decay;
      

      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  auto dist = (ip->distributions);
	  // set the proton momenta in ps
	  // -> needed to define the boost to lab-frame
	  //    in FillDistributions()
	  ps->P1() = 1.0/x[0] * (ps->p1());
	  ps->P2() = 1.0/x[1] * (ps->p2());
	  	      
	  double vwgt  = ip->cmp_v_weight();
	  if (flags & F_EVAL_B_QCDxQCD)
	    { // only QCD LO: fill in the 0-th hostogram
	      ps->FillDistributions(*dist,H_LO_QCD,res_b*vwgt,mScale);
// #ifdef WITH_T_SPIN
// 	      FV const& k1 = ps->k1();
// 	      FV const& k2 = ps->k2();
// 	      (*dist)[9]->FillOne(H_LO_QCD,obs_M12(k1,k2)*mScale,res_b*vwgt);
// #endif
	    }
	  else
	    { // PHI + INT LO: fill in the first histogram
	      ps->FillDistributions(*dist,H_LO_PHI,res_b*vwgt,mScale);
// #ifdef WITH_T_SPIN
// 	      FV const& k1 = ps->k1();
// 	      FV const& k2 = ps->k2();
// 	      c_double const& D1 = hm->GetBoson(0)->GetPropagator();
// 	      c_double const& D2 = hm->GetBoson(1)->GetPropagator();
// 	      c_double D1D2 = D1*std::conj(D2);
// 	      (*dist)[7]->FillOne(H_LO_PHI,obs_M12(k1,k2)*mScale,res_b*vwgt*D1D2.real());
// 	      (*dist)[9]->FillOne(H_LO_PHI,obs_M12(k1,k2)*mScale,res_b*vwgt);
// #endif
	    }
	  if (EVAL_V(flags)) ps->FillDistributions(*dist,H_NLO_PHI_V,res_v*vwgt,mScale);//res_v*
	}
      ////////////////////////////////////////////////////////////////////////////
      return (res_b+res_v);
    }
  else
    {
      return 0.0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
double Integrand_2_2_pdf_ID(double* x, size_t dim, void* arg)
{
  using namespace Constants;

#ifdef WITH_T_SPIN
  if (unlikely(dim!=14))
#else
    if (unlikely(dim!=4))
#endif
      {
	ERROR("wrong integral dimension!");
      }

  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_2*         ps    = dynamic_cast<PS_2_2*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   MUF2   = hm->MUF2();
  double const&   mt2    = hm->mt2();
  double const&   mScale = hm->Scale(); 
  double const&   mScale2= hm->Scale2();
  
  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;
  // static const double cf_qg = PREF_QG*CONV_GeV2i_pbarn/mScale2;
  
  // hadronic c.m.e.
  double& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];

  
  double res_d = 0.0;
  double res_dx_gg = 0.0;
  double res_dx_qg = 0.0;
  double res_dx_gq = 0.0;
  if (ps->set(sqrt(s_part),x[2],0.0))
    {
      // gluon PDFs
      // scale MU must be given in [GeV] !!! 
      double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*s_part);
      
      // finite part of the integrated dipoles: delta terms and distribution end-point terms
      res_d = Eval_ID_GG(*ps,x[3],*hm,flags)*jac_flux_1*cf_gg*f1_g*f2_g;
#ifdef DEBUG
      CHECKNA(res_d);
#endif
      
      // DISTRIBUTIONS: x = 1 contribution ////////////////////////////////////////
      if (ip->collect_dist)
	{
	  auto dist = (ip->distributions);
	  // set the proton momenta in ps
	  // -> needed to define the boost to lab-frame
	  //    in FillDistributions()
	  ps->P1() = 1.0/x[0] * (ps->p1());
	  ps->P2() = 1.0/x[1] * (ps->p2());

	  double vwgt  = ip->cmp_v_weight();
	  ps->FillDistributions(*dist,H_NLO_PHI_ID,(res_d)*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////   

      
      // boosted phase space (boost is done in Eval_ID_X)
      static PS_2_2 ps_x1(mt2,mt2,"x*p1 p2 -> k1' k2'  [static in Integrand_2_2_pdf_ID]");
      static PS_2_2 ps_x2(mt2,mt2,"p1 x*p2 -> k1' k2'  [static in Integrand_2_2_pdf_ID]");
  
      // finite part of the integrated dipoles: x-dependent terms of P and K operators
      double sx = x[3]*s_part;
      if ( ps_x1.set(sqrt(sx),x[2],0.0) )
	{
	  if (Eval_ID_X(ps_x1,x[3],*hm,flags,res_dx_gg,res_dx_qg))
	    {
	      // modified version of flux and phase space density
	      double jac_flux_x = TwoPi*ps_x1.get_wgt()/(2.0*sx);
	      // GG contribution
	      res_dx_gg *= jac_flux_x*cf_gg*f1_g*f2_g;
	      // QG/GQ contribution
	      if (flags & F_EVAL_D_QG_CONT)
		{
		  // need quark PDFs
		  static auto quark_pids     = { 5,  4,  3,  2,  1};	  
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
		  res_dx_qg *= jac_flux_x*cf_gg*(qg_pdf+gq_pdf);
		  // the following is not necessary, Eval_ID_X_QG() is invariant under p1<->p2
		  // (this is mainly because the gg->tt Born MEs are symmetric under p1<->p2,
		  //  while the IF/FI colour correlated Born MEs are asymmetric)
		  // ps_x1.swap_initial_state();
		  // res_dx_gq += Eval_ID_X_QG(ps_x1,*hm,flags)*jac_flux_x*cf_gg*gq_pdf;
		  // ps_x1.swap_initial_state();
		}
	      // DISTRIBUTIONS: x neq 1 contribution //////////////////////////////////////
	      if (ip->collect_dist)
		{
		  static LT dip_boost;
		  auto dist = (ip->distributions);
		  double vwgt  = ip->cmp_v_weight();
		  ps_x1.P1() = (ps->P1());
		  ps_x1.P2() = (ps->P2());
		  ps_x1.set_x(x[3]);
		  ps_x2 = ps_x1;
	      
		  dip_boost.set_boost_z(x[3],0);
		  dip_boost.apply(ps_x1.k1());
		  dip_boost.apply(ps_x1.k2());
	      
		  dip_boost.set_boost_z(x[3],1);
		  dip_boost.apply(ps_x2.k1());
		  dip_boost.apply(ps_x2.k2());
	      
		  ps_x1.FillDistributions(*dist,H_NLO_PHI_ID,0.5*(res_dx_gg+res_dx_qg)*vwgt,mScale);
		  ps_x2.FillDistributions(*dist,H_NLO_PHI_ID,0.5*(res_dx_gg+res_dx_qg)*vwgt,mScale);
		}
	      ////////////////////////////////////////////////////////////////////////////  
	      
	    }
	}
    }
  return  (res_d+res_dx_gg+res_dx_qg+res_dx_gq);
}








double Integrand_2_3_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  
  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_3*         ps    = dynamic_cast<PS_2_3*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   MUF2   = hm->MUF2();
  double const&   mScale = hm->Scale(); 
  double const&   mScale2= hm->Scale2();

  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;
  
  // hadronic c.m.e.
  double const& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];

  // make a technical cut on the angle between p1/p2 and p3 
  if ( ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5]))
    {
      if (fabs(x[2])  > Cuts::COLL_CUT) return 0.0;
      if (ps->p3()[0] < Cuts::SOFT_CUT) return 0.0;
      
// #ifdef WITH_T_SPIN
//       // // top anti-top 4-vectors given in the tt z.m.f.
//       // FV const& k1 = ps->k1();
//       // FV const& k2 = ps->k2();
//       FV& S1 = ps->s1();
//       FV& S2 = ps->s2();
//       S1 = ps->s1_r();
//       S2 = ps->s2_r();
      
//       // boost from k1 restframe to tt z.m.f.
//       static LT boost_k1_RF;
//       boost_k1_RF.set_boost(k1,1);
//       boost_k1_RF.apply(S1);
//       // boost from k2 restframe to tt z.m.f.
//       static LT boost_k2_RF;
//       boost_k2_RF.set_boost(k2,1);
//       boost_k2_RF.apply(S2);
// #endif

      // the PDF wants dimensionful quantities in units of GeV (here, MUF is in units of mt!) 
      double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_gg = TwoPi*ps->get_wgt()/(2.0*s_part)*cf_gg*f1_g*f2_g;

      double res_r_gg = Eval_R_GG(*ps,*hm,flags)*jac_flux_gg;

      // for the unint. dipoles we need here already the VEGAS weight
      double vwgt = 0.0;
      std::vector<HistArray*>* dist = nullptr;
      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  dist = (ip->distributions);
	  ps->P1() = 1.0/x[0] * (ps->p1());
	  ps->P2() = 1.0/x[1] * (ps->p2());
	  vwgt  = ip->cmp_v_weight();
	  ps->FillDistributions(*dist,H_NLO_PHI_R,res_r_gg*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////
      
      double res_d_gg = Eval_UID_GG (*ps,*hm,flags,-vwgt*jac_flux_gg,dist)*jac_flux_gg;

      return (res_r_gg-res_d_gg);
    }
  else
    {
      return 0.0;
    }
}



double Integrand_2_3_qg_qq_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  //using namespace RunParameters;
  
  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_3*         ps    = dynamic_cast<PS_2_3*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   MUF2   = hm->MUF2();
  double const&   mScale = hm->Scale(); 
  double const&   mScale2= hm->Scale2();

  // spin/color average for gg initial state, conversion GeV^-2 -> picobarn
  static const double cf_gg = PREF_GG*CONV_GeV2i_pbarn/mScale2;  
  // spin/color average for qg/gq initial state, conversion GeV^-2 -> picobarn
  static const double cf_qg = PREF_QG*CONV_GeV2i_pbarn/mScale2;
  // spin/color average for qq initial state, conversion GeV^-2 -> picobarn
  static const double cf_qq = PREF_QQ*CONV_GeV2i_pbarn/mScale2;
  
  // hadronic c.m.e.
  double s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];
  
  // make a technical cut on the angle between p1/p2 and p3 
  if ( ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5]))
    {
      // technical cuts on soft/collinear phase space
      if (fabs(x[2])  > Cuts::COLL_CUT) return 0.0;
      if (ps->p3()[0] < Cuts::SOFT_CUT) return 0.0;
      
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);
      // need VEGAS weight explicitly for the uint. dipoles later
      double vwgt = 0.0;
      std::vector<HistArray*>* dist = nullptr;
      /////////////////////////// PDF //////////////////////////////////
      double qq_pdf = 0.0;
      double qg_pdf = 0.0;
      double gq_pdf = 0.0;
      double f1_g = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2_g = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
      // get the quark/antiquark pdfs
      static auto quark_pids     = { 5,  4,  3,  2,  1};	  
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
      
      double res_r_qq = Eval_R_QQ(*ps,*hm,flags)*cf_qq*jac_flux*qq_pdf;
      double res_r_qg = Eval_R_QG(*ps,*hm,flags)*cf_qg*jac_flux*qg_pdf;

      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  dist = (ip->distributions);
	  if (unlikely(dist==nullptr))
	    {
	      ERROR("You want distributions? Then pass something more than a nullptr.");
	    }
	  ps->P1() = 1.0/x[0] * (ps->p1());
	  ps->P2() = 1.0/x[1] * (ps->p2());
	  vwgt  = ip->cmp_v_weight();
	  // fill qq and qg weights
	  ps->FillDistributions(*dist,H_NLO_PHI_R,(res_r_qq+res_r_qg)*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////
      double res_d_qg = Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_gg*jac_flux*qg_pdf,dist)*cf_gg*jac_flux*qg_pdf;
  
      ////////////////////////////////////////////////////////////////////////////
      //// crossed phase space ///////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////
      ps->swap_initial_state();
      double res_r_gq = Eval_R_QG(*ps,*hm,flags)*cf_qg*jac_flux*gq_pdf;
      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  ps->FillDistributions(*dist,H_NLO_PHI_R,res_r_gq*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////
      double res_d_gq = Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_gg*jac_flux*gq_pdf,dist)*cf_gg*jac_flux*gq_pdf;
      ps->swap_initial_state();    
      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////

      // if (obs_PT12(ps->k1(),ps->k2())*mScale < 0.1)
      // 	{
      // 	  std::cout << "\n=================================================================";
      // 	  PRINT(obs_PT12(ps->k1(),ps->k2())*mScale);
      // 	  PRINT(res_r_qq);
      // 	  PRINT(res_r_qg);
      // 	  PRINT(res_d_qg);
      // 	  PRINT(res_r_gq);
      // 	  PRINT(res_d_gq);
      // 	  std::cout << "=================================================================\n";
      // 	}
      return res_r_qq + (res_r_qg-res_d_qg) + (res_r_gq-res_d_gq);
    }
  else
    {
      return 0.0;
    }
}
