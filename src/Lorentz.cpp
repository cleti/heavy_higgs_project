
#include "../inc/Lorentz.h"




////// class FV     ///////////////////////////////////////////////////////////////////////////////////////
FV& FV::operator=(std::initializer_list<double> L)
{
  if (L.size() == 4)
    {
      v = L;
    }
  else
    {
      ERROR("expecting initializer list of length 4");
    }
  return *this;
}
FV& FV::operator=(FV const& other)
{
  //std::cout << std::endl << " FV: assignment (copy)" << std::endl;
  if (this!=&other)
    {
      v = other.v;
    }
  return *this;
}
FV& FV::operator=(FV&& other) noexcept
{
  //std::cout << std::endl << " FV: assignment (move)" << std::endl;
  if (this!=&other)
    {
      v = other.v;
      // v.swap(other.v);
    }
  return *this;
}

FV& FV::operator+=(FV const& other)
{
  for (int i=0;i<4;++i)
    {
      v[i]+=other.v[i];
    }
  return *this;
}
FV& FV::operator-=(FV const& other)
{
  for (int i=0;i<4;++i)
    {
      v[i]-=other.v[i];
    }
  return *this;
}
FV& FV::operator*=(double const& a)
{
  for (int i=0;i<4;++i)
    {
      v[i]*=a;
    }
  return *this;
}
FV& FV::operator/=(double const& a)
{
  if (a==0)
    {
      ERROR("division by 0");
    }
  for (int i=0;i<4;++i)
    {
      v[i]/=a;
    }
  return *this;
}
void FV::swap(FV& other)
{
  v.swap(other.v);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
FV operator*(FV const& v, double const& a)
{
  //std::cout << std::endl << " FV: binary operator * (copy). " << std::endl;
  FV V = v;
  V*=a;
  return V;
}
FV operator*(double const& a, FV const& v)
{
  //std::cout << std::endl << " FV: binary operator * (copy). " << std::endl;
  FV V = v;
  V*=a;
  return V;
}
FV&& operator*(double const& a, FV&& v)
{
  //std::cout << std::endl << " FV: binary operator * (copy). " << std::endl;
  return std::move(v*=a);
}
FV&& operator*(FV&& v, double const& a)
{
  //std::cout << std::endl << " FV: binary operator * (copy). " << std::endl;
  return std::move(v*=a);
}
////////////////////////////////////////////////////////////////////////////////////////////////


FV operator+(FV const& v1, FV const& v2)
{
  //std::cout << std::endl << " FV: binary operator + (copy). " << std::endl;
  FV V = v1;
  V+=v2;
  return V;
}
FV&& operator+(FV const& v1, FV&& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator + (move 1). " << std::endl;
  return std::move(v2+=v1);
}
FV&& operator+(FV&& v1,FV const& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator + (move 2). " << std::endl;
  return std::move(v1+=v2);
}
FV&& operator+(FV&& v1,FV&& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator + (move 3). " << std::endl;
  return std::move(v1+=v2);
}
////////////////////////////////////////////////////////////////////////////////////////////////
FV operator-(FV const& v1,FV const& v2)
{
  //std::cout << std::endl << " FV: binary operator - (copy). " << std::endl;
  FV V = v1;
  V-=v2;
  return V;
}
FV&& operator-(FV const& v1, FV&& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator - (move 1). " << std::endl;
  return std::move(operator-(std::move(v2-=v1)));
}
FV&& operator-(FV&& v1,FV const& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator - (move 2). " << std::endl;
  return std::move(v1-=v2);
}
FV&& operator-(FV&& v1,FV&& v2) noexcept
{
  //std::cout << std::endl << " FV: binary operator - (move 3). " << std::endl;
  return std::move(v1-=v2);
}
////////////////////////////////////////////////////////////////////////////////////////////////
FV operator-(FV const& v1)
{
  //std::cout << std::endl << " FV: unary operator - (copy). " << std::endl;
  FV V = v1;
  V*=(-1);
  return V;
}
FV&& operator-(FV&& v1) noexcept
{
  //std::cout << std::endl << " FV: unary operator - (move). " << std::endl;
  return std::move(v1*=(-1));
}
////////////////////////////////////////////////////////////////////////////////////////////////





////// class LT     ///////////////////////////////////////////////////////////////////////////////////////
const double LT::G[4] = {1.0,-1.0,-1.0,-1.0};

LT::LT()
{}
LT::~LT()
{}

void LT::transpose()
{
  for (int i=0;i<4;++i)
    {
      for (int j=i;j<4;++j)
	{
	  std::swap(M[i][j],M[j][i]);
	}
    }
}

void LT::invert()
{
  for (int i=1;i<4;++i)
    {
      if (M[0][i] != 0.0)
	{
	  M[0][i] *= (-1);
	  M[i][0] *= (-1);
	}
    }
}

int LT::set_FF(FV const& p1, FV const& p2)
{
  double p1_p1 = sp(p1,p1);
  double p2_p2 = sp(p2,p2);
  double p1_p2 = sp(p1,p2);


  double t2 = p1_p1-p2_p2;
  double t3 = p1_p1+p2_p2;
  double t4 = p1_p2*p1_p2;
  double t5 = p1_p1*p2_p2;
  double t6 = 1.0/(2.0*t5);
  double t7 = sqrt(t4-t5);

  double sinhx = t6*(-t2*p1_p2+t3*t7);
  double coshx = t6*( t3*p1_p2-t2*t7);

  double A = sinhx/t7;
  double B = (coshx-1.0)/(t7*t7);


#ifdef DEBUG_LT
  std::cout << std::endl;
#endif
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{
	  M[i][j] = A*(p1[i]*p2[j]-p1[j]*p2[i]) + B*(p1_p2*(p1[i]*p2[j]+p1[j]*p2[i])-p1_p1*p2[i]*p2[j]-p2_p2*p1[i]*p1[j]);
	  if (i==j) M[i][j] += G[i];
#ifdef DEBUG_LT
	  std::cout << std::setw(8) << M[i][j] << " | ";
#endif
	}
#ifdef DEBUG_LT
      std::cout << std::endl;
#endif
    }
#ifdef DEBUG_LT
  std::cout << std::endl;
#endif

#ifdef DEBUG_LT
  std::cout << std::endl;
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{

	  double test = 0.0;
	  for (int k=0;k<4;++k)
	    {
	      test += M[i][k]*M[j][k]*G[k];
	    }

	  std::cout << std::setw(8) << test << " | ";
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;
#endif

  return 1;
}

void LT::apply(FV& v)
{
  double temp[4] = {v[0],v[1],v[2],v[3]};
  for (int i=0;i<4;++i)
    {
      v[i] = 0.0;
      for (int j=0;j<4;++j)
	{
	  v[i] += M[i][j]*temp[j];
	}
    }
}
void LT::apply_G(FV& v)
{
  double temp[4] = {v[0],v[1],v[2],v[3]};
  for (int i=0;i<4;++i)
    {
      v[i] = 0.0;
      for (int j=0;j<4;++j)
	{
	  v[i] += M[i][j]*temp[j]*G[j];
	}
    }
}

int LT::set_wigner(FV const& P1, FV const& P2)
{
  double m1 = sqrt(MSQ(P1));
  double g1 = P1[0]/m1;
  if ( std::isnan(g1) || std::isinf(g1) )
    {
      WARNING("gamma_1 is inf/NaN, maybe m_1 is too small!");
      PRINT(m1);
      return 0;
    }
  double beta1[4] = {1.0,P1[1]/P1[0],P1[2]/P1[0],P1[3]/P1[0]};

  double m2 = sqrt(MSQ(P2));
  double g2 = P2[0]/m2;
  if ( std::isnan(g2) || std::isinf(g2) )
    {
      WARNING("gamma_2 is inf/NaN, maybe m_2 is too small!");
      PRINT(m2);
      return 0;
    }
  double beta2[4] = {1.0,P2[1]/P2[0],P2[2]/P2[0],P2[3]/P2[0]};

  double b1b2 = beta1[1]*beta2[1]+beta1[2]*beta2[2]+beta1[3]*beta2[3];
  double g = 1.0+g1*g2*(1.0+b1b2);


  double g1_1m = 1.0-g1;
  double g1_1p = 1.0+g1;
  double g2_1m = 1.0-g2;
  double g2_1p = 1.0+g2;
  double g1g2 = g1*g2;

  M[0][0] = 1.0;
  for (int i=1;i<4;++i)
    {
      M[i][0] = 0.0;
      M[0][i] = 0.0;
      for (int j=1;j<4;++j)
	{

	  M[i][j] = g1g2/(g)*( beta1[i]*beta1[j]*g1/g2*(g2_1m)/(g1_1p) + beta2[i]*beta2[j]*g2/g1*(g1_1m)/(g2_1p) - beta1[i]*beta2[j] + beta2[i]*beta1[j]*(1.0+2.0*g1g2/(g1_1p*g2_1p)*b1b2) );
	  if (i==j) M[i][j]+=1.0;
	}
    }

  return 1;
}


// transformation of final state momenta for contruction of unintegrated initial-initial dipoles
int LT::set_II(FV const& K, FV const& Kb)
{

  double t1 = sp(K,K);
  double t2 = sp(Kb,Kb);
  double t3 = 2.0*sp(K,Kb);
  double A  = 2.0/(t1+t2+t3);
  double B  = 2.0/t1;


#ifdef DEBUG_LT
  std::cout << std::endl;
#endif
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{
	  M[i][j] = -A*(K[i]+Kb[i])*(K[j]+Kb[j])+B*Kb[i]*K[j];	  
	  if (i==j) M[i][j] += G[i];
#ifdef DEBUG_LT
	  std::cout << std::setw(8) << M[i][j] << " | ";
#endif
	}
#ifdef DEBUG_LT
      std::cout << std::endl;
#endif
    }
#ifdef DEBUG_LT
  std::cout << std::endl;
#endif

#ifdef DEBUG_LT
  std::cout << std::endl;
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{

	  double test = 0.0;
	  for (int k=0;k<4;++k)
	    {
	      test += M[i][k]*M[j][k]*G[k];
	    }

	  std::cout << std::setw(8) << test << " | ";
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;
#endif

  return 1;
}



// boost to the restframe of the given 4-vector
// setting INV=true inverts the boost, i.e. vectors given in the restframe of P
// will be transformed to the current frame
int LT::set_boost(FV const& P,bool INV)
{
  double m = sqrt(MSQ(P));
  double gamma = P[0]/m;
  if ( std::isnan(gamma) || std::isinf(gamma) )
    {
      WARNING("gamma is inf/NaN, maybe m is too small!");
      PRINT(m);
      return 0;
    }
  // here one should check beta[i] > 1e-8;
  double beta[4] = {1.0,P[1]/P[0],P[2]/P[0],P[3]/P[0]};
  double beta2 = beta[1]*beta[1]+beta[2]*beta[2]+beta[3]*beta[3];
  if ( std::isnan(beta2) || std::isinf(beta2) )
    {
      WARNING("beta is inf/NaN, maybe P[0] is too small!");
      PRINT(P[0]);
      return 0;
    }


  M[0][0] = gamma;
  for (int i=1;i<4;++i)
    {
      M[0][i] = -beta[i]*gamma;
      if (INV)
	{
	  M[0][i] *= (-1);
	}
      M[i][0] = M[0][i];
      for (int j=i;j<4;++j)
	{
	  // if (i==0)
	  //   {
	  //     M[i][j] = beta[j-1]*gamma;
	  //     M[j][i] = M[i][j];
	  //   }
	  // else
	  //   {
	  if (beta2>1e-15) M[i][j] = (gamma-1.0)*beta[i]*beta[j]/beta2;
	  else M[i][j] = 0.0;
	      if ( std::isnan(M[i][j])  )
		{
		  WARNING("matrix component is NaN!");
		  std::cout << std::setprecision(16);
		  std::cout << " component [" << i << "," << j << "]"  << std::endl;
		  std::cout << " b = [" << beta[0] << " " << beta[1] << " " << beta[2] << "], b^2 = " << beta2 << ", g = " << gamma << std::endl;
		  SLEEP(3);
		  return 0;
		}
		if (i==j)
		  {
		    M[i][j] += 1.0;
		  }
		else
		  {
		    M[j][i] = M[i][j];
		  }
		
	    // }
	}
    }


#ifdef DEBUG_LT
  print();

  std::cout << std::endl;
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{

	  double test = 0.0;
	  for (int k=0;k<4;++k)
	    {
	      test += M[i][k]*M[j][k]*G[k];
	    }

	  std::cout << std::setw(8) << test << " | ";
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;
#endif
  return 1;
}

// transforms 4-vectors from the z.m.f. of x*p1, p2
// to the z.m.f. of p1,p2 (provieded p1,p2 in z-direction)
// if INV=true: tranformation from p1, x*p2 z.m.f. to p1,p2 z.m.f.
int LT::set_boost_z(double const& x,bool INV)
{
  // x = sqrt((1+beta)/(1-beta))
  // -> beta = (1-x^2)/(1+x^2)
  if (x<=0.0 || x>=1.0)
    {
      WARNING("x out of range! provide 0<x<1!");
      return 0;
    }
 
  double beta = (x-1.0)/(1.0+x);
  double gamma= (1.0+x)/(2.0*sqrt(x));

  M[0][0] = gamma;
  M[0][1] = 0.0;
  M[0][2] = 0.0;
  M[0][3] = beta*gamma;
  if (INV) M[0][3] *= (-1);

  M[1][1] = 1.0;
  M[1][2] = 0.0;
  M[1][3] = 0.0;  

  M[2][2] = 1.0;
  M[2][3] = 0.0;  

  M[3][3] = gamma;

  
  for (int i=0;i<4;++i)
    {
      for (int j=i+1;j<4;++j)
	{
	  M[j][i] = M[i][j];
	}
    }
  return 1;
}


void LT::print()
{
  std::cout << std::endl;
  for (int i=0;i<4;++i)
    {
      for (int j=0;j<4;++j)
	{
	  std::cout << std::setw(20) << M[i][j] << " | ";
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
