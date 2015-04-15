

#ifndef GLOBAL_H
#define GLOBAL_H

#include <complex>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <iostream>
#include <streambuf>

#include "Makros.h"

typedef std::complex<double> c_double;
typedef unsigned long ulong;

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Constants {
  /* misc. constants */
  extern const double Pi;
  extern const double gE;
  extern const double Pi2;
  extern const double Pi3;
  extern const double gE2;
  extern const double TwoPi;
  extern const double FourPi;
  extern const double TwoPi2;
  /* conversion factors */
  extern const double CONV_MeV_fm;
  extern const double CONV_GeV2i_mbarn;
  extern const double CONV_GeV2i_pbarn;
  /* first two coefficients of the expansion of (4 pi)^eps/Gamma(1-eps) */
  extern const double C_eps1;
  extern const double C_eps2;
  /* colour factors */
  extern const double CA;
  extern const double TF;
  extern const double CF;
  extern const double CFCA2;
  extern const double Nf;
  extern const double CA2;
  extern const double TF2;
  extern const double CF2;
  extern const double beta0;
  /* initial gg colour and spin average */
  extern const double PREF_GG;
  extern const double PREF_QQ;
  extern const double PREF_QG;
  /* top mass */
  extern const double MT;

  extern const double BR_TT_LL;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
namespace Cuts {
  extern double COLL_CUT;
  extern double SOFT_CUT;
  extern double IDIP_X_CUT;
}
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
// teebuf and teestream
class teebuf: public std::streambuf
{
 public:
  // Construct a streambuf which tees output to both input
  // streambufs.
 teebuf(std::streambuf * sb1, std::streambuf * sb2) :
  sb1(sb1),
  sb2(sb2)
 { }
 private:
  // This tee buffer has no buffer. So every character "overflows"
  // and can be put directly into the teed buffers.
  virtual int overflow(int c)
  {
    if (c == EOF)
      {
	return !EOF;
      }
    else
      {
	int const r1 = sb1->sputc(c);
	int const r2 = sb2->sputc(c);
	return (r1 == EOF || r2 == EOF) ? EOF : c;
      }
  }
    
  // Sync both teed buffers.
  virtual int sync()
  {
    int const r1 = sb1->pubsync();
    int const r2 = sb2->pubsync();
    return (r1 == 0 && r2 == 0) ? 0 : -1;
  }   
 private:
  std::streambuf * sb1;
  std::streambuf * sb2;
};

class teestream : public std::ostream
{
 public:
  // Construct an ostream which tees output to the supplied
  // ostreams.
 teestream(std::ostream & o1, std::ostream & o2) :
  std::ostream(&tbuf),
  tbuf(o1.rdbuf(), o2.rdbuf())
 { }      
 private:
  teebuf tbuf;
};
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#endif
