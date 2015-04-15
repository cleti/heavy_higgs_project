

#include "../inc/Integrands_pp_ttX.h"


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
  //using namespace RunParameters;

  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_2*         ps    = dynamic_cast<PS_2_2*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   mScale = hm->Scale();  
  double const&   mScale2= hm->Scale2();
   
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
      double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*rs_part*rs_part);
      
      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      //// no initial gluon helicity average ! only the sum !
      double cf = 4.0*PREF_GG*CONV_GeV2i_pbarn/mScale2;

      // born matrix elements (GG)
      double res_b = Eval_B(*ps,*hm,flags,0)*jac_flux_1*cf;
      
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

  if (unlikely(dim!=3))
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
  
  // hadronic c.m.e.
  double& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];
  
  if (ps->set(sqrt(s_part),x[2],0.0))
    {
      
#ifdef WITH_T_SPIN
      if (ip->tDecay)
	{
	  /////////////////////////////////////////////////////////////////////////////////////
	  // evaluate decay phase spaces and matrix elements
	  /////////////////////////////////////////////////////////////////////////////////////
	  double tt_decay = 1.0;
	  static PS_2_3* ps_t1 = static_cast<PS_2_3*>(ps->get_child(0));
	  static PS_2_3* ps_t2 = static_cast<PS_2_3*>(ps->get_child(1));

	  if (unlikely(ps_t1 == nullptr)) {ERROR("'ps_t1': nullptr received");}
	  if (unlikely(ps_t2 == nullptr)) {ERROR("'ps_t2': nullptr received");}

	  // in the end this sould be the top/antitop spin vectors in the p1+p2 z.m.f. 
	  FV& S1 = ps->s1();
	  FV& S2 = ps->s2();
		  
	  if ( ps_t1->set(mt2,x[4] ,x[5] ,x[6] ,x[7] ,x[8] ,0.0,0.0,0.0) &&
	       ps_t2->set(mt2,x[9] ,x[10],x[11],x[12],x[13],0.0,0.0,0.0)    )
	    {
	      ///////////////// decay of top quark ////////////////////////////////////////////
	      tt_decay *= (ps_t1->get_wgt() * Eval_t_blnu(*ps_t1,S1) / (GammaT));
	      // spin vector in top restframe
	      /////////////////////////////////////////////////////////////////////////////////

	      ///////////////// decay of anti-top quark ///////////////////////////////////////
	      tt_decay *= (ps_t2->get_wgt() * Eval_t_blnu(*ps_t2,S2) / (GammaT));
	      S2 *= (-1);// spin vector in antitop restframe
	      /////////////////////////////////////////////////////////////////////////////////
	      // correcting the phase space factors of 2->2 * ( 2->3 * 2->3 )
	      // to match that of 2->6 phase space
	      tt_decay *= mt2/Pi2; 
	    }
	  else
	    {
	      return 0.0;
	    }
	}
      ////////////////////////////////////////////////////////////////////
      // boost top spin and decay vectors to parton z.m.f. ///////////////
      ////////////////////////////////////////////////////////////////////
      // top anti-top 4-vectors given in the tt z.m.f.
      FV const& k1 = ps->k1();
      FV const& k2 = ps->k2();
      FV const& S1_r = ps->s1_r();
      FV const& S2_r = ps->s2_r();
      FV& S1 = ps->s1();
      FV& S2 = ps->s2();
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
      ////////////////////////////////////////////////////////////////////
#endif
      
      double res_b = 0.0;
      double res_v = 0.0;
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*s_part);

      // gluon PDFs
      // scale MU must be given in [GeV] !!! 
      double f1 = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2 = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];

      // include spin/color average of initial gluons and PFDs, convert to units of picobarn
      double cf = PREF_GG*f1*f2*CONV_GeV2i_pbarn/mScale2;
      
      // born matrix elements (GG)
      if (EVAL_B(flags)) res_b = Eval_B(*ps,*hm,flags,0)*jac_flux_1*cf;

      // born matrix elements (QQ)
      if (flags & F_EVAL_B_QCDxQCD)
        {
          static auto quark_pids     = { 5,  4,  3,  2,  1};	  
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
          res_b += Eval_B_QQ(*ps,*hm)*jac_flux_1*PREF_QQ*(qq_pdf)*CONV_GeV2i_pbarn/mScale2; 
        }
	  
      // finite part of the virtual corrections (GG)
      if (EVAL_V(flags)) res_v = Eval_V(*ps,*hm,flags)*jac_flux_1*cf;
      

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
	    }
	  else
	    { // PHI + INT LO: fill in the first histogram
	      ps->FillDistributions(*dist,H_LO_PHI,res_b*vwgt,mScale);
	    }
	  if (EVAL_V(flags)) ps->FillDistributions(*dist,H_NLO_V,res_v*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////
      return  (res_b+res_v);
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
  //using namespace RunParameters;

  if (unlikely(dim!=4))
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
  
  // hadronic c.m.e.
  double& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double  s_part = s_hadr*x[0]*x[1];

  
  double res_d = 0.0;
  double res_dx= 0.0;

  // gluon PDFs
  // scale MU must be given in [GeV] !!! 
  double f1 = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
  double f2 = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
  // include spin/color average of initial gluons and PFDs, convert to units of picobarn
  double cf = PREF_GG*f1*f2*CONV_GeV2i_pbarn/mScale2*BR_TT_LL;

  if (ps->set(sqrt(s_part),x[2],0.0))
    {
      ps->set_x(x[3]);
      
#ifdef WITH_T_SPIN
      // top anti-top 4-vectors given in the tt z.m.f.
      FV const& k1 = ps->k1();
      FV const& k2 = ps->k2();
      FV const& S1_r = ps->s1_r();
      FV const& S2_r = ps->s2_r();
      FV& S1 = ps->s1();
      FV& S2 = ps->s2();
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
#endif
 
      // flux factor 1/(2*s_part), factor 2 pi from trivial integration over phi
      double jac_flux_1 = TwoPi*ps->get_wgt()/(2.0*s_part);

      // finite part of the integrated dipoles: delta terms and distribution end-point terms
      res_d = Eval_ID(*ps,*hm,flags)*jac_flux_1*cf;
#ifdef DEBUG
      CHECKNA(res_d);
#endif



      // boosted phase space (boost is done in Eval_ID_X)
      static PS_2_2 ps_x(mt2,mt2,"p1 p2 -> k1 k2 (boosted) [static in Integrand_2_2_pdf_ID]");
  
      // finite part of the integrated dipoles: x-dependent terms of P and K operators
      double sx = x[3]*s_part;

      if ( ps_x.set(sqrt(sx),x[2],0.0) )
	{
	  ps_x.set_x(x[3]);
      
	  // modified version of flux and phase space density
	  res_dx = Eval_ID_X(ps_x,*hm,flags);
      
	  // the ps_x object is modified in Eval_ID_X -> use the pahsespace weight only after applying this function !!!
	  double jac_flux_x = TwoPi*ps_x.get_wgt()/(2.0*sx);
#ifdef DEBUG
	  CHECKNA(jac_flux_x);
#endif
      
	  res_dx *= jac_flux_x*cf;
#ifdef DEBUG
	  CHECKNA(res_dx);
#endif
	}
  
      // DISTRIBUTIONS: x neq 1 contribution //////////////////////////////////////
      if (ip->collect_dist)
	{
	  auto dist = (ip->distributions);
	  // set the proton momenta in ps
	  // -> needed to define the boost to lab-frame
	  //    in FillDistributions()
	  ps->P1() = 1.0/x[0] * (ps->p1());
	  ps->P2() = 1.0/x[1] * (ps->p2());
	  ps_x.P1() = (ps->P1());
	  ps_x.P2() = (ps->P2());

	  double vwgt  = ip->cmp_v_weight();
	  ps->FillDistributions(*dist,H_NLO_ID,(res_d)*vwgt,mScale);
	  ps_x.FillDistributions(*dist,H_NLO_ID,res_dx*vwgt,mScale);
	}
      ////////////////////////////////////////////////////////////////////////////   
    }


  
  return  (res_d+res_dx);
}








double Integrand_2_3_pdf(double* x, size_t dim, void* arg)
{
  using namespace Constants;
  //using namespace RunParameters;
  
  integrand_par*  ip    = static_cast<integrand_par*> (arg);
  PS_2_3*         ps    = dynamic_cast<PS_2_3*> (ip->ps);
  const ulong&    flags = ip->eval_flags;
  HiggsModel*     hm    = ip->higgs_model;
  
  double const&   MUF2   = hm->MUF2();
  double const&   mScale2= hm->Scale2();
  
  // hadronic c.m.e.
  double const& s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];

  // make a technical cut on the angle between p1/p2 and p3 
  if ( ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5]))
    {
      if (fabs(x[2])>Cuts::COLL_CUT) return 0.0;
      FV const& p3 = ps->p3();
      if (p3[0]<Cuts::SOFT_CUT) return 0.0;
      
#ifdef WITH_T_SPIN
      // // top anti-top 4-vectors given in the tt z.m.f.
      // FV const& k1 = ps->k1();
      // FV const& k2 = ps->k2();
      FV& S1 = ps->s1();
      FV& S2 = ps->s2();
      S1 = ps->s1_r();
      S2 = ps->s2_r();
      
      // boost from k1 restframe to tt z.m.f.
      static LT boost_k1_RF;
      boost_k1_RF.set_boost(k1,1);
      boost_k1_RF.apply(S1);
      // boost from k2 restframe to tt z.m.f.
      static LT boost_k2_RF;
      boost_k2_RF.set_boost(k2,1);
      boost_k2_RF.apply(S2);


#endif

      
      double jac_flux = TwoPi*ps->get_wgt()/(2.0*s_part);

      // the PDF wants dimensionful quantities in units of GeV (here, MUF is in units of mt!) 
      double f1 = ip->pdf->xfxQ2(21, x[0], MUF2*mScale2) / x[0];
      double f2 = ip->pdf->xfxQ2(21, x[1], MUF2*mScale2) / x[1];
      
      double cf_gg = PREF_GG*f1*f2*CONV_GeV2i_pbarn/mScale2*BR_TT_LL;
      double res_r_gg = Eval_R_GG(*ps,*hm,flags)*cf_gg*jac_flux;
      
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
	  ps->FillDistributions(*dist,H_NLO_R,res_r_gg*vwgt);
	}
      ////////////////////////////////////////////////////////////////////////////
      
      double res_d_gg = Eval_UID_GG (*ps,*hm,flags,-vwgt*cf_gg*jac_flux,dist)*cf_gg*jac_flux;

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
  double const&   mScale2= hm->Scale2();
  
  // hadronic c.m.e.
  double s_hadr = ip->s_hadr;
  // partonic c.m.e.
  double s_part = s_hadr*x[0]*x[1];

  // make a technical cut on the angle between p1/p2 and p3 
  if ( ps->set( sqrt(s_part),x[2],0.0,x[3],x[4],x[5]))// && fabs(x[2])<COLL_CUT )
    {
      if (fabs(x[2])>Cuts::COLL_CUT) return 0.0;
      FV const& p3 = ps->p3();
      if (p3[0] < Cuts::SOFT_CUT) return 0.0;

      //// SPIN STUFF 
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
      
      double cf_qq = PREF_QQ*qq_pdf*CONV_GeV2i_pbarn/mScale2;
      double cf_qg = PREF_QG*qg_pdf*CONV_GeV2i_pbarn/mScale2;
      double cf_gq = PREF_QG*gq_pdf*CONV_GeV2i_pbarn/mScale2;
      
      double res_r_qq = Eval_R_QQ(*ps,*hm,flags)*cf_qq*jac_flux;
      double res_r_qg = Eval_R_QG(*ps,*hm,flags)*cf_qg*jac_flux;
     
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
	  ps->FillDistributions(*dist,H_NLO_R,(res_r_qq+res_r_qg)*vwgt);
	}
      ////////////////////////////////////////////////////////////////////////////
      double res_d_qg = Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_qg*jac_flux,dist)*cf_qg*jac_flux;

      
      //////////////////////////////////////////////////////////////////////////
      //// crossed phase space ///////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////
      ps->swap();
      double res_r_gq = Eval_R_QG(*ps,*hm,flags)*cf_gq*jac_flux;
      // DISTRIBUTIONS ///////////////////////////////////////////////////////////
      if (ip->collect_dist)
	{
	  // dist = (ip->distributions);
	  // if (unlikely(dist==nullptr))
	  //   {
	  //     ERROR("You want distributions? Then pass something more than a nullptr.");
	  //   }
	  // ps->P1() = 1.0/x[0] * (ps->p1());
	  // ps->P2() = 1.0/x[1] * (ps->p2());
	  // vwgt  = ip->cmp_v_weight();
	  // fill gq weight
	  ps->FillDistributions(*dist,H_NLO_R,res_r_gq*vwgt);
	}
      ////////////////////////////////////////////////////////////////////////////
      double res_d_gq = Eval_UID_QG (*ps,*hm,flags,-vwgt*cf_gq*jac_flux,dist)*cf_gq*jac_flux;
      ps->swap();    
      ////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////
      
      return res_r_qq + (res_r_qg-res_d_qg) + (res_r_gq-res_d_gq);
    }
  else
    {
      return 0.0;
    }
}
