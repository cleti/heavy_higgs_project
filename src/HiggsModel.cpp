#include "../inc/HiggsModel.h"




//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
HiggsBoson::HiggsBoson(
		       double const& M,
		       double const& G,
		       double const& Vh,
		       double const& at,
		       double const& bt,
		       double const& ab,
		       double const& bb) :
  d_M(M),
  d_M2(M*M),
  d_G(G),
  d_G2(G*G),
  d_At(at/Vh),
  d_Bt(bt/Vh),
  d_Ab(ab/Vh),
  d_Bb(bb/Vh),
  d_FH_eff( d_At/12.0),
  d_FA_eff(bt?(-d_Bt/16.0):0),
  d_S_FF(0),
  d_FH_full(0),
  d_FA_full(0),
  d_S_Den(0),
  d_DenSq(0),
  d_Den(0)
{
  if (Vh<=0.0)
    {
      std::cout << std::endl << " In " << __FUNCTION__ << ", " << __LINE__ << std::endl;
      std::cout << " ERROR: value of Vh needs to be non-zero!" << std::endl;
      d_At = 0;
      d_Bt = 0;
      d_Ab = 0;
      d_Bb = 0;
      d_FH_eff = 0;
      d_FA_eff = 0;
    }
}


void HiggsBoson::SetFormFactors(
				double const& S,
				double const& mt2,
				double const& mb2)
{
  // we only need the finite part of the itnegrals
  const int I = 0;
  const double zero = 0.0;
  // as long as mt>0, mb>0 the integrals are finite and do not depend on MUR
  const double MUR2 = 0.0;
  
  if (S!=d_S_FF)
    {
      d_S_FF = S;
      // reset values
      d_FH_full = 0.0;
      d_FA_full = 0.0;
	
      // // top contribution
      // if (d_At != 0.0 || d_Bt != 0.0)
      // 	{
      c_double I3t = qli3_(&S,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I);
      d_FH_full += mt2 * d_At * ( (S-4.0*mt2) * I3t - 2.0 );
      d_FA_full += mt2 * d_Bt * S * I3t;
      // }
      // mb = 0 indicates that b-quark effects on the effective gg-phi coupling should be negledted
      if (mb2>0.0)
	{
	  // // bottom contribution
	  // if (d_Ab != 0.0 || d_Bb != 0.0)
	  //   {
	  c_double I3b = qli3_(&S,&zero,&zero,&mb2,&mb2,&mb2,&MUR2,&I);
	  d_FH_full += mb2 * d_Ab * ( (S-4.0*mb2) * I3b - 2.0 );
	  d_FA_full += mb2 * d_Bb * S * I3b;
	  // }
	}
      d_FH_full /= -(4.0*S);
      d_FA_full /=  (8.0*S);

#ifdef DEBUG
      if (S<=0.0)
	{
	  ERROR("S<=0!");
	}   
      if (std::isnan(RE(d_FH_full)) || std::isinf(RE(d_FH_full)) || std::isnan(IM(d_FH_full)) || std::isinf(IM(d_FH_full)))
	{
	  PRINT(qli3_(&S,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I));
	  PRINT(qli3_(&S,&zero,&zero,&mb2,&mb2,&mb2,&MUR2,&I));	  
	  ERROR("d_FH_full is inf/nan!");
	}
      if (std::isnan(RE(d_FA_full)) || std::isinf(RE(d_FA_full)) || std::isnan(IM(d_FA_full)) || std::isinf(IM(d_FA_full)))
	{
	  PRINT(qli3_(&S,&zero,&zero,&mt2,&mt2,&mt2,&MUR2,&I));
	  PRINT(qli3_(&S,&zero,&zero,&mb2,&mb2,&mb2,&MUR2,&I));	   
	  ERROR("d_FA_full is inf/nan!");
	}
#endif
    }
}

void HiggsBoson::SetPropagator(double const& S, bool rescale)
{
  if (S!=d_S_Den)
    {
      d_S_Den = S;
      d_DenSq = 1.0 / (pow(S-d_M2,2)+pow(d_M*d_G,2));
      d_Den   = c_double((S-d_M2)*d_DenSq,-d_M*d_G*d_DenSq);

#ifdef DEBUG
      if (S<=0.0)
	{
	  ERROR("S<=0!");
	}        
      if (std::isnan(d_DenSq) || std::isinf(d_DenSq))
	{
	  ERROR("d_DenSq is inf/nan!");
	}
      if (std::isnan(RE(d_Den)) || std::isinf(RE(d_Den)) || std::isnan(IM(d_Den)) || std::isinf(IM(d_Den)))
	{
	  ERROR("d_Den is inf/nan!");
	}
#endif
      // rescale effective Higgs-gluon vertex by ratio full vertex/eff. vertex
      // here it is assumed that each scalar propagator is connected to one ggH vertex
      // assumes that d_FH_full and d_FA_full have been set up correctly
      if (rescale)
      	{
      	  //c_double K1 = ((d_FH_full)+4.0*(d_FA_full))/((d_FH_eff)+4.0*(d_FA_eff));
      	  double   K2 = (std::norm(d_FH_full)+4.0*std::norm(d_FA_full))/(std::norm(d_FH_eff)+4.0*std::norm(d_FA_eff));
	  std::cout << std::endl << K2 << std::endl;
      	  // this affects the interference terms QCD x PHI_i as well as PHI_i x PHI_j
      	  d_Den   *= std::sqrt(K2);
      	  // this affects only the squared terms PHI_i x PHI_i
      	  d_DenSq *= K2;
      	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void HiggsPrefactors::Reset()
{
  At_fH_re = 0.0;		
  At_fA_re = 0.0;		
  Bt_fH_re = 0.0;
  Bt_fA_re = 0.0;
  At_fH_im = 0.0;
  At_fA_im = 0.0;
  Bt_fH_im = 0.0;
  Bt_fA_im = 0.0;
  At2_fH2_De = 0.0;
  At2_fA2_De = 0.0;
  Bt2_fH2_De = 0.0;
  Bt2_fA2_De = 0.0;
  At_Bt_fH2_De = 0.0;
  At_Bt_fA2_De = 0.0;
  At_Bt_fH2_DeIM = 0.0;
  At_Bt_fA2_DeIM = 0.0;
}
void HiggsPrefactors::Print(std::ostream& ost)
{
  SPRINT(ost,At_fH_re);		
  SPRINT(ost,At_fA_re);		
  SPRINT(ost,Bt_fH_re);
  SPRINT(ost,Bt_fA_re);
  SPRINT(ost,At_fH_im);
  SPRINT(ost,At_fA_im);
  SPRINT(ost,Bt_fH_im);
  SPRINT(ost,Bt_fA_im);
  SPRINT(ost,At2_fH2_De);
  SPRINT(ost,At2_fA2_De);
  SPRINT(ost,Bt2_fH2_De);
  SPRINT(ost,Bt2_fA2_De);
  SPRINT(ost,At_Bt_fH2_De);
  SPRINT(ost,At_Bt_fA2_De);
  SPRINT(ost,At_Bt_fH2_DeIM);
  SPRINT(ost,At_Bt_fA2_DeIM);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void AmplitudePrefactors::Print(std::ostream& ost)
{
/* born */
SPRINT(ost,PREF_B_PHIxPHI);
SPRINT(ost,PREF_B_PHIxQCD);
SPRINT(ost,PREF_B_QCDxQCD);
SPRINT(ost,PREF_B_QCDxQCD_CF);
SPRINT(ost,PREF_B_QCDxQCD_CA);
SPRINT(ost,PREF_B_QCDxQCD_CFCA);
  /* virtual */
SPRINT(ost,PREF_V_PHIxQCD);
SPRINT(ost,PREF_V_PHIxQCD_CF);
SPRINT(ost,PREF_V_PHIxQCD_CA);
SPRINT(ost,PREF_V_PHIxQCD_CFCA2);
SPRINT(ost,PREF_V_PHIxQCD_Nf);
SPRINT(ost,PREF_V_PHIxQCD_CT);
SPRINT(ost,PREF_V_PHIxPHI);
SPRINT(ost,PREF_V_PHIxPHI_CA);
SPRINT(ost,PREF_V_PHIxPHI_CF);
  /* real */
SPRINT(ost,PREF_R_PHIxQCD);
SPRINT(ost,PREF_R_PHIxQCD_CF);
SPRINT(ost,PREF_R_PHIxQCD_CA);
SPRINT(ost,PREF_R_PHIxQCD_CFCA2);
SPRINT(ost,PREF_R_PHIxPHI);
SPRINT(ost,PREF_R_PHIxPHI_CA);
SPRINT(ost,PREF_R_PHIxPHI_CF);
  /* UID */
SPRINT(ost,PREF_UID_TF);
SPRINT(ost,PREF_UID_CA);
SPRINT(ost,PREF_UID_CF);  
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
HiggsModel::HiggsModel(std::string const & name):
  d_Name(name),
  d_AlphaS(0),
  d_AlphaS2(0),
  d_MUR(0),
  d_MUR2(0),
  d_MUF(0),
  d_MUF2(0),
  d_mt(0),
  d_mt2(0),
  d_mb(0),
  d_mb2(0),
  d_VH(1),
  d_Scale(1),
  d_Scale2(1),
  d_useK(false)
{

}
HiggsModel::~HiggsModel()
{

}


void HiggsModel::SetHiggsPrefactors(double const& S, bool EFF)
{
  d_HPref.Reset();
  // iterate over d_Bosons and add the individual contributions
  auto last  = std::end(d_Bosons);
  // PRINT(std::sqrt(S));
  for (auto phi_i = std::begin(d_Bosons); phi_i!=last; ++phi_i)
    {
      if (!EFF || d_useK) (*phi_i)->SetFormFactors(S,d_mt2,d_mb2);
      // multiply prop. with K-factor only if effective vertices are used
      // in this case we also need to compute the full form factors
      // when the full vertices are used in the amplitudes do nothing
      (*phi_i)->SetPropagator(S,d_useK && EFF);
      c_double const& FHi = (*phi_i)->GetFH(EFF);
      c_double const& FAi = (*phi_i)->GetFA(EFF);
      c_double const& Di  = (*phi_i)->GetPropagator();
      double FH2i  = std::norm(FHi);
      double FA2i  = std::norm(FAi);
      double const& di  = (*phi_i)->GetPropagatorSq();
      double const& Ati = (*phi_i)->At();
      double const& Bti = (*phi_i)->Bt();
      double At2i = std::pow(Ati,2);
      double Bt2i = std::pow(Bti,2);

      // PRINT(Ati);
      // PRINT(Bti);
      // PRINT(Di*FHi*(-4.0*S));
      // PRINT(Di*FAi*( 8.0*S));
      // PRINT(Di);

      
      // prefactors in PHIxQCD amplitudes
      // these are just the sum of individual phi contributions
      d_HPref.At_fH_re += Ati*(FHi*Di).real();
      d_HPref.At_fA_re += Ati*(FAi*Di).real();
      d_HPref.Bt_fH_re += Bti*(FHi*Di).real();
      d_HPref.Bt_fA_re += Bti*(FAi*Di).real();

      d_HPref.At_fH_im += Ati*(FHi*Di).imag();
      d_HPref.At_fA_im += Ati*(FAi*Di).imag();
      d_HPref.Bt_fH_im += Bti*(FHi*Di).imag();
      d_HPref.Bt_fA_im += Bti*(FAi*Di).imag();

      // prefactors in PHIxPHI amplitudes
      // first the sum of the individual phi_i^2 contributions
      d_HPref.At2_fH2_De += At2i*FH2i*di;
      d_HPref.At2_fA2_De += At2i*FA2i*di;
      d_HPref.Bt2_fH2_De += Bt2i*FH2i*di;
      d_HPref.Bt2_fA2_De += Bt2i*FA2i*di;
	  
#ifdef WITH_T_SPIN
      // only needed for polarized amplitudes:
      // PHIxPHI terms: add at*bt interference from Phi1^2 contribution
      d_HPref.At_Bt_fH2_De += Ati*Bti*FH2i*di;
      d_HPref.At_Bt_fA2_De += Ati*Bti*FA2i*di;   
#endif

      // second loop to get the phi^2 contributions
      for (auto phi_j = phi_i+1; phi_j!=last; ++phi_j)
	{
	  if (!EFF || d_useK) (*phi_j)->SetFormFactors(S,d_mt2,d_mb2);
	  // multiply prop. with K-factor only if effective vertices are used
	  // in this case we also need to compute the full form factors
	  // when the full vertices are used in the amplitudes do nothing
	  (*phi_j)->SetPropagator(S,d_useK && EFF);

	  c_double const& FHj = (*phi_j)->GetFH(EFF);
	  c_double const& FAj = (*phi_j)->GetFA(EFF);
	  c_double const& Dj  = (*phi_j)->GetPropagator();
	  double const& Atj = (*phi_j)->At();
	  double const& Btj = (*phi_j)->Bt();
      
	  c_double Dij  = Di*std::conj(Dj);
	  c_double FHij_Dij = FHi*std::conj(FHj)*Dij;
	  c_double FAij_Dij = FAi*std::conj(FAj)*Dij;
	  
	  // prefactors in PHIxPHI amplitudes
	  // here we also have to sum all the phi_i * phi_j interferences    
	  d_HPref.At2_fH2_De   += 2.0*Ati*Atj*(FHij_Dij).real();
	  d_HPref.At2_fA2_De   += 2.0*Ati*Atj*(FAij_Dij).real();
	  d_HPref.Bt2_fH2_De   += 2.0*Bti*Btj*(FHij_Dij).real();
	  d_HPref.Bt2_fA2_De   += 2.0*Bti*Btj*(FAij_Dij).real();
#ifdef WITH_T_SPIN
	  // only needed for polarized amplitudes:
	  /*
	    What is the purpose of this function?
	    In the PHIxPHI amplitudes with 2 (or more) Higgs bosons we encounter interferences
	    
	      RE[D1*(D2^*) * A1*(A2^*)] =
	          RE[D1*(D2^*)] * RE[A1*(A2^*)] - IM[D1*(D2^*)] * IM[A1*(A2^*)] = 

		  RE[D1*(D2^*)] * (at1*at2*RE[fH1*(fH2^*)]*RE[A(a,H)*(A(a,H)^*)] + ... ) -
		  RE[D1*(D2^*)] * (at1*at2*IM[fH1*(fH2^*)]*IM[A(a,H)*(A(a,H)^*)] + ... ) -
	          IM[D1*(D2^*)] * (at1*at2*RE[fH1*(fH2^*)]*IM[A(a,H)*(A(a,H)^*)] + ... ) -
		  IM[D1*(D2^*)] * (at1*at2*IM[fH1*(fH2^*)]*RE[A(a,H)*(A(a,H)^*)] + ... ) =

		  ( RE[D1*(D2^*)]*at1*at2*RE[fH1*(fH2^*)] - 
		    IM[D1*(D2^*)]*at1*at2*IM[fH1*(fH2^*)] ) * RE[A(a,H)*(A(a,H)^*)] + ...  -
		  ( RE[D1*(D2^*)]*at1*at2*IM[fH1*(fH2^*)] + 
		    IM[D1*(D2^*)]*at1*at2*RE[fH1*(fH2^*)] ) * IM[A(a,H)*(A(a,H)^*)] - ...  =		  
		    at1*at2*RE[fH1*(fH2^*)*D1*(D2^*)] * RE[A(a,H)*(A(a,H)^*)] + ...  
		  - at1*at2*IM[fH1*(fH2^*)*D1*(D2^*)] * IM[A(a,H)*(A(a,H)^*)] - ...  ,
		    
	    where Ai are the parts of the respective gg(qq/qg) -> phi (g/q) -> tt(g/q)
	    amplitudes without the propagator.
	    Extracting all Higgs-specific coupling factors, the first term gives
	    
	      RE[D1*(D2^*)] * RE[A1*(A2^*)] =
	        RE[D1*(D2^*)] * (at1*at2*RE[fH1*(fH2^*)]*RE[A(a,H)*(A(a,H)^*)] + ... ) -
		RE[D1*(D2^*)] * (at1*at2*IM[fH1*(fH2^*)]*IM[A(a,H)*(A(a,H)^*)] + ... )
		
	    ( ... = all further combinations of ati,bti,fHi,fAi).
	    The amplitudes RE[A(a,H)*(A(a,H)^*)] are independent of all specific Higgs couplings.
	    They have already been computed for the case of a single Higgs boson
	      
	      RE[D1*(D1^*) * A1*(A1^*)] =
	        |D1|^2 * (at1^2*at1*|fH1|^2*RE[A(a,H)*(A(a,H)^*)] + ... ),

	    we just have to pick up the coefficients of D1*at1*|fH1|^2 ... .
	    Thus by replacing the prefactors

	      |D1|^2 * (at1^2*at1*|fH1|^2 ->
	      
	        |D1|^2 * (at1^2*at1*|fH1|^2 +
		|D2|^2 * (at2^2*at2*|fH2|^2 +
		RE[D1*(D2^*)] * at1*at2*RE[fH1*(fH2^*)]*RE[A(a,H)*(A(a,H)^*)],

	      ...
		
	    in these amplitudes we get all necesssary terms
	    except those proportional to IM[A1*(A2^*).
	    For the imaginary parts we have a similar decomposition
	    
	      IM[D1*(D2^*)] * IM[A1*(A2^*)] =
	        IM[D1*(D2^*)] * (at1*at2*RE[fH1*(fH2^*)]*IM[A(a,H)*(A(a,H)^*)] + ... ) +
		IM[D1*(D2^*)] * (at1*at2*IM[fH1*(fH2^*)]*RE[A(a,H)*(A(a,H)^*)] + ... ),

	    but in our case only the terms IM[A(a,H/A)*(A(b,H/A)^*)] are nonzero
	    and only in the case when we do not sum over top spins.
		

		
	  */
	  // PHIxPHI terms: at Phi1 * bt  Phi2 interferences
	  d_HPref.At_Bt_fH2_De += (Ati*Btj+Atj*Bti)*(FHij_Dij).real();
	  d_HPref.At_Bt_fA2_De += (Ati*Btj+Atj*Bti)*(FAij_Dij).real();	
	  // the Phi1 At * Phi2 Bt interference is prop. to the imaginary parts
	  d_HPref.At_Bt_fH2_DeIM -= (Ati*Btj-Atj*Bti)*(FHij_Dij).imag();
	  d_HPref.At_Bt_fA2_DeIM -= (Ati*Btj-Atj*Bti)*(FAij_Dij).imag();
#endif
	}
    }
  // EXIT(1);
}


void HiggsModel::SetAmpPrefactors()
{
    using namespace Constants;
    /* born */
    d_APref.PREF_B_PHIxPHI      = CF*CA2*d_AlphaS2/Pi2;
    d_APref.PREF_B_PHIxQCD      = CF*CA*d_AlphaS2;
    d_APref.PREF_B_QCDxQCD      = CF*CA*d_AlphaS2*Pi2;
    d_APref.PREF_B_QCDxQCD_CF   = CF2*CA*d_AlphaS2*Pi2;
    d_APref.PREF_B_QCDxQCD_CA   = CF*CA2*d_AlphaS2*Pi2;
    d_APref.PREF_B_QCDxQCD_CFCA = CF*CA*CFCA2*d_AlphaS2*Pi2;
    /* virtual */
    d_APref.PREF_V_PHIxQCD_CF    = d_APref.PREF_B_PHIxQCD*d_AlphaS/Pi*CF;
    d_APref.PREF_V_PHIxQCD_CA    = d_APref.PREF_B_PHIxQCD*d_AlphaS/Pi*CA;
    d_APref.PREF_V_PHIxQCD_CFCA2 = d_APref.PREF_B_PHIxQCD*d_AlphaS/Pi*CFCA2;
    d_APref.PREF_V_PHIxQCD_Nf    = d_APref.PREF_B_PHIxQCD*d_AlphaS/Pi*Nf;
    d_APref.PREF_V_PHIxQCD_CT    = d_APref.PREF_B_PHIxQCD*d_AlphaS/Pi;
    d_APref.PREF_V_PHIxQCD       = d_APref.PREF_B_PHIxQCD*d_AlphaS;
    d_APref.PREF_V_PHIxPHI       = d_APref.PREF_B_PHIxPHI*d_AlphaS/Pi;
    d_APref.PREF_V_PHIxPHI_CA    = d_APref.PREF_B_PHIxPHI*d_AlphaS/Pi*CA;
    d_APref.PREF_V_PHIxPHI_CF    = d_APref.PREF_B_PHIxPHI*d_AlphaS/Pi*CF;
    /* real */
    d_APref.PREF_R_PHIxQCD       = d_APref.PREF_V_PHIxQCD_CT*Pi2;
    d_APref.PREF_R_PHIxQCD_CF    = d_APref.PREF_V_PHIxQCD_CF*Pi2;
    d_APref.PREF_R_PHIxQCD_CA    = d_APref.PREF_V_PHIxQCD_CA*Pi2;
    d_APref.PREF_R_PHIxQCD_CFCA2 = d_APref.PREF_V_PHIxQCD_CFCA2*Pi2;
    d_APref.PREF_R_PHIxPHI       = d_APref.PREF_V_PHIxPHI*Pi2;
    d_APref.PREF_R_PHIxPHI_CA    = d_APref.PREF_V_PHIxPHI_CA*Pi2;
    d_APref.PREF_R_PHIxPHI_CF    = d_APref.PREF_V_PHIxPHI_CF*Pi2;
    /* UID */
    d_APref.PREF_UID_TF =  4.0*Pi*d_AlphaS;
    d_APref.PREF_UID_CA = 16.0*Pi*d_AlphaS*CA;
    d_APref.PREF_UID_CF =  8.0*Pi*d_AlphaS*CF;
    /* Needed in  EVAL_V  */
    d_APref.MUR2   = d_MUR2;
    d_APref.AlphaS = d_AlphaS;
}

// void SetScale(double const& val)
// {
//   // // first rescale all values with the previous scale

// }


void HiggsModel::AddBoson(
			  double const& M,
			  double const& G,
			  double const& a_t,
			  double const& b_t,
			  double const& a_b,
			  double const& b_b)
{
  if (d_VH>0.0)
    {
      d_Bosons.push_back(std::make_shared<HiggsBoson>(M/d_Scale,
						      G/d_Scale,
						      d_VH,// already normalized to d_Scale
						      a_t,b_t,a_b,b_b));
    }
  else
    {
      ERROR("call SetVH(val) to set up the combined VEV");
    }
}

void HiggsModel::AddScalar(
			   double const& M,
			   double const& G,
			   double const& a_t,
			   double const& a_b)
{
  if (d_VH>0.0)
    {
      d_Bosons.push_back(std::make_shared<HiggsBoson>(M/d_Scale,
						      G/d_Scale,
						      d_VH,// already normalized to d_Scale
						      a_t,0.0,a_b,0.0));
    }
  else
    {
      ERROR("call SetVH(val) to set up the combined VEV");
    }
}

void HiggsModel::AddPseudoscalar(
		     double const& M,
		     double const& G,
		     double const& b_t,
		     double const& b_b)
{
  if (d_VH>0.0)
    {
      d_Bosons.push_back(std::make_shared<HiggsBoson>(M/d_Scale,
						      G/d_Scale,
						      d_VH,// already normalized to d_Scale
						      0.0,b_t,0.0,b_b));
    }
  else
    {
      ERROR("call SetVH(val) to set up the combined VEV");
    }
}
void HiggsModel::PopBoson()
{
  d_Bosons.pop_back();
}



void HiggsModel::Print(std::ostream& ost, double const& mScale) const
{
  ost << std::endl << std::setfill('=') << std::setw(100) << "" << std::setfill(' ');
  ost << std::endl << " Model configuration [" << d_Name << "]: " << std::endl;
  SPRINT(ost,d_AlphaS);
  SPRINT(ost,d_MUR*mScale);  
  SPRINT(ost,d_MUF*mScale);
  SPRINT(ost,d_mt*mScale);
  SPRINT(ost,d_mb*mScale);
  SPRINT(ost,d_VH*mScale);
  SPRINT(ost,d_Scale);
  
  ost << std::endl << " Higgs bosons:" << std::endl << std::endl;
  if (d_Bosons.size()>0)
    {
      ost << std::setw(5) <<  "#" << std::setw(15) <<  "M"  << std::setw(15) <<  "Gamma" << std::setw(15) <<  "a_t" << std::setw(15) <<  "b_t" << std::setw(15) <<  "a_b" << std::setw(15) <<  "b_b" << std::endl;
      ost << std::setfill('-') << std::setw(100) << "" << std::setfill(' ') << std::endl;
      int I = 1;
      for (auto phi:d_Bosons)
	{
	  ost << std::setw(5) << I << std::setw(15) <<  phi->M()*mScale << std::setw(15) << phi->G()*mScale << std::setw(15) << phi->At()*d_VH << std::setw(15) <<  phi->Bt()*d_VH << std::setw(15) << phi->Ab()*d_VH << std::setw(15) <<  phi->Bb()*d_VH << std::endl;
	  ++I;
	}
      // d_HPref.Print(ost);
    }
  else
    {
      ost << " no Higgs bosons added so far. " << std::endl;
    }
  ost << std::setfill('=') << std::setw(100) << "" << std::setfill(' ') << std::endl;
  // d_APref.Print(ost);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
