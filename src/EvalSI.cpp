

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <complex>

#include <gsl/gsl_sf_dilog.h>

#include "../inc/Global.h"
#include "../inc/Makros.h"


#ifdef __cplusplus
extern "C" {
#endif
  extern void qlinit_();
  extern std::complex<double> qli1_(double*,double*,int*);
  extern std::complex<double> qli2_(double*,double*,double*,double*,int*);
  extern std::complex<double> qli3_(double*,double*,double*,double*,double*,double*,double*,int*);
  extern std::complex<double> qli4_(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*);
#ifdef __cplusplus
}
#endif

int main()
{
  qlinit_(); 
  std::cout << std::endl;
  std::cout.flush();

  int P = 10;
  int W = 15;

  double x1,x2,x3,x4;
  double s1,s2;
  double m0,m1,m2,m3;
  double mu2=1.0;
  int    i;
  std::complex<double> res(0.0,0.0);


  {
    mu2 = 2.5;

    
    m0 = 1.0;

    std::cout << std::endl << "ql_i1(mt^2) = ";
    i = -2;
    res = qli1_(&m0,&mu2,&i);
    c_double I1_2 = res;
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli1_(&m0,&mu2,&i);
    c_double I1_1 = res;
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli1_(&m0,&mu2,&i);
    c_double I1_0 = res;
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;


    double s = 15.0;
    x1  = 1.0;
    x2  = s;
    x3  = 1.0;

    m1  = 0.0;
    m2  = 1.0;

    std::cout << std::endl << "ql_i3(m^2,s,m^2,0,m^2,m^2) = ";
    i = -2;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3a_2 = res;
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3a_1 = res;
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3a_0 = res;
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;

    x1  = 0.0;
    x2  = s;
    x3  = 0.0;

    m1  = 0.0;
    m2  = 0.0;

    std::cout << std::endl << "ql_i3(0,s,0,0,0,0) = ";
    i = -2;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3b_2 = res;
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3b_1 = res;
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli3_(&x1,&x2,&x3,&m1,&m2,&m2,&mu2,&i);
    c_double I3b_0 = res;
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;


    
    double beta2 = (1.0-4.0/s);
    double beta  = sqrt(beta2);
    {
      c_double S0 = -0.5*(1.0+beta2)*s*I3a_0 - (1.0-beta2)*s*I3a_1 - I1_0;
      c_double P0 = -0.5*(1.0+beta2)*s*I3a_0 - I1_0;

      c_double H0 = 2.0*s*I3b_0 - 11.0/3.0;
      c_double A0 = H0 - 4.0 + 11.0/3.0;

      double ItS  = 2.0*CF*S0.imag()-CA*H0.imag();
      double ItP  = 2.0*CF*P0.imag()-CA*H0.imag();
      double Rt1S = 2.0*CF*S0.real()-CA*H0.real();
      double Rt2S = 2.0*CF*S0.real()-CA*A0.real();
      double Rt1P = 2.0*CF*P0.real()-CA*A0.real();
      double Rt2P = 2.0*CF*P0.real()-CA*H0.real();

      PRINT(S0);
      PRINT(P0);
      PRINT(H0);
      PRINT(A0);
	      
      PRINT(ItS);
      PRINT(ItP);
      PRINT(Rt1S);
      PRINT(Rt2S);
      PRINT(Rt1P);
      PRINT(Rt2P);

      std::cout << std::endl;
    }

    {
      double x  = (1.0-beta)/(1.0+beta);
      double x1 =  1.0-x;
      double c1 =  0.5*(1.0+beta2)/beta;
      double c2 =      (1.0-beta2)/beta;
      
      c_double S0 = c_double( c1*( 2.0/3.0*Pi2 - 0.5*pow(log(x),2) + 2.0*(log(x1)*log(x)+dilog(x)) ) - c2*log(x) - 1.0 ,
			      -Pi*c1*( log(x) - 2.0*log(x1) + c2/c1 )  );
      c_double P0 = c_double( S0.real() + c2*log(x) , S0.imag() + Pi*c2 );

      double RH = pow(log(s/mu2),2) - Pi2 - 11.0/3.0;
      double RA = RH - 4.0 + 11.0/3.0;
      double IG = -TwoPi*log(s/mu2);
      
      c_double EPS0 = c_double(-2.0*CF*(c1*log(x)+1.0)*log(mu2) , -2.0*CF*Pi*c1*log(mu2) );

      double ItS  = 2.0*CF*S0.imag()-CA*IG+EPS0.imag();
      double ItP  = 2.0*CF*P0.imag()-CA*IG+EPS0.imag();
      double Rt1S = 2.0*CF*S0.real()-CA*RH+EPS0.real();
      double Rt2S = 2.0*CF*S0.real()-CA*RA+EPS0.real();
      double Rt1P = 2.0*CF*P0.real()-CA*RA+EPS0.real();
      double Rt2P = 2.0*CF*P0.real()-CA*RH+EPS0.real();

      PRINT(S0+EPS0/(2.0*CF));
      PRINT(P0+EPS0/(2.0*CF));

      PRINT(RH);
      PRINT(RA);
      PRINT(IG);

      PRINT(EPS0);
      
      PRINT(ItS);
      PRINT(ItP);
      PRINT(Rt1S);
      PRINT(Rt2S);
      PRINT(Rt1P);
      PRINT(Rt2P);
      std::cout << std::endl;
    }

    
    exit(1);









    
    mH2 = pow(500.0,2)/mScale2;
    
    x1 = mH2;
    x2 = 0.0;
    x3 = 0.0;
    
    m2 = mt2;
    mu2 = mH2;

    std::cout << std::endl << "ql_i3(s,0,0,m^2,m^2,m^2) = ";
    i = -2;
    res = qli3_(&x1,&x2,&x3,&m2,&m2,&m2,&mu2,&i);
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli3_(&x1,&x2,&x3,&m2,&m2,&m2,&mu2,&i);
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli3_(&x1,&x2,&x3,&m2,&m2,&m2,&mu2,&i);
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;

    c_double test = c_double(0.0,0.0);
    if (4.0*m2>mH2)
      {
	test = c_double(-2.0/x1*pow(std::asin(1.0/sqrt(4.0*m2/x1)),2),0.0);
      }
    else
      {
	c_double beta = sqrt(1.0-4.0*m2/mH2);
	test = 1.0/(2.0*x1)*pow(std::log((beta+1.0)/(beta-1.0)),2);
      }
    std:: cout << std::endl << " check = " << test << std::endl;
    
    exit(1);
    
    x3 = 0.0;
    x1 = 0.0;
    x2 = 12.234;
             
    m0 = 0.0;
    m1 = 0.0;
    m2 = 0.0;
    mu2 = 9.23;

    std::cout << std::endl << "2*s*ql_i3(0,0,s,0,0,0) = ";
    i = -2;
    res = qli3_(&x1,&x2,&x3,&m0,&m1,&m2,&mu2,&i);
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << 2.0*x2*res.real() << " + i " << std::setprecision(P) << std::setw(W) << 2.0*x2*res.imag() << " ]";
    i = -1;
    res = qli3_(&x1,&x2,&x3,&m0,&m1,&m2,&mu2,&i);
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << 2.0*x2*res.real() << " + i " << std::setprecision(P) << std::setw(W) << 2.0*x2*res.imag() << " ]";
    i = 0;
    res = qli3_(&x1,&x2,&x3,&m0,&m1,&m2,&mu2,&i);
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << 2.0*x2*res.real() << " + i " << std::setprecision(P) << std::setw(W) << 2.0*x2*res.imag() << " ]" << std::endl;

    std::cout << std::endl << "ln(s/mu2)^2-pi^2 = " << pow(log(x2/mu2),2)-Pi2 << std::endl;
    std::cout << std::endl << "-2*pi*ln(s/mu2)  = " << -TwoPi*log(x2/mu2) << std::endl;

    exit(1);

    x1 = 5.0;
    m0 = 0.0;
    m1 = 0.0;

    std::cout << std::endl << "ql_i2(s,0,0) = ";
    i = -2;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;


    x1 = 1.0;
    m0 = 1.0;
    m1 = 0.0;

    std::cout << std::endl << "ql_i2(mt^2,mt^2,0) = ";
    i = -2;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "        eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = -1;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "      + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
    i = 0;
    res = qli2_(&x1,&m0,&m1,&mu2,&i);
    std::cout << std::endl << "      + eps^0  * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;


    double MT2 = 1.0;
    double S12 = 5.0;
    double T11 = -1.0;
    // double T12 = -2.0;

    x1 = MT2;
    x2 = S12;
    x3 = 0.0;
    x4 = T11;
    s1 = MT2;
    s2 = 0.0;

    m0 = 0.0;
    m1 = MT2;
    m2 = MT2;
    m3 = MT2;

      std::cout << std::endl << "ql_i4 = ";
      i = -2;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << "  eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
      i = -1;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << " + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
      i = 0;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << " + eps^0 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;


    x1 = 0.0;
    x2 = 0.0;
    x3 = MT2;
    x4 = MT2;
    s1 = S12;
    s2 = T11;

    m0 = MT2;
    m1 = MT2;
    m2 = MT2;
    m3 = 0.0;

      std::cout << std::endl << "ql_i4 = ";
      i = -2;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << "  eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
      i = -1;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << " + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
      i = 0;
      res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
      std::cout << " + eps^0 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;



  }


  // for (int j=0;j<1;j++)
  //   {
  //     x1 = x1_vals[j];
  //     x2 = x2_vals[j];
  //     x3 = x3_vals[j];
  //     x4 = x4_vals[j];
  //     s1 = s1_vals[j];
  //     s2 = s2_vals[j];
  //     m0 = m0_vals[j];
  //     m1 = m1_vals[j];
  //     m2 = m2_vals[j];
  //     m3 = m3_vals[j];

  //     std::cout << std::endl << "{ "<< x1 << ", " << x2 << ", " << x3 << ", " << x4 <<  "; "<< s1 << ", " << s2 << "; " << m0 << ", " << m1 << ", " << m2 << ", " << m3  << " }" << std::endl;

  //     std::cout << std::endl << "ql_i4 = ";
  //     i = -2;
  //     res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
  //     std::cout << "  eps^-2 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
  //     i = -1;
  //     res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
  //     std::cout << " + eps^-1 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]";
  //     i = 0;
  //     res = qli4_(&x1,&x2,&x3,&x4,&s1,&s2,&m0,&m1,&m2,&m3,&mu2,&i);
  //     std::cout << " + eps^0 * [" << std::setprecision(P) << std::setw(W) << res.real() << " + i " << std::setprecision(P) << std::setw(W) << res.imag() << " ]" << std::endl;

  //   }



  std::cout << std::endl;
  exit(1);
}
