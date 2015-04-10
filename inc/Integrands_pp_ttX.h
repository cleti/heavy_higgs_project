

#ifndef INTEGRANDS_PPTT_H
#define INTEGRANDS_PPTT_H

#include "Global.h"
#include "Functions_Shared.h"

#include "Functions_pp_ttX_V.h"
#include "Functions_pp_ttX_ID.h"
#include "Functions_pp_ttX_R.h"
#include "Functions_pp_ttX_UID.h"


double Integrand_poly2(double* x, size_t dim, void* arg);

// integrand without PDFs to reproduce Bernies s_part plot
double Integrand_2_2(double* x, size_t dim, void* arg);

double Integrand_2_2_pdf_BV(double* x, size_t dim, void* arg);

double Integrand_2_2_pdf_ID(double* x, size_t dim, void* arg);

double Integrand_2_3_pdf(double* x, size_t dim, void* arg);

double Integrand_2_3_qg_qq_pdf(double* x, size_t dim, void* arg);

#endif
