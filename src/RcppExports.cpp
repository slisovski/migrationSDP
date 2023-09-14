// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Predation
double Predation(const int& time, const int& site, const int& x, double u, double f);
RcppExport SEXP _migrationSDP_Predation(SEXP timeSEXP, SEXP siteSEXP, SEXP xSEXP, SEXP uSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const int& >::type site(siteSEXP);
    Rcpp::traits::input_parameter< const int& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(Predation(time, site, x, u, f));
    return rcpp_result_gen;
END_RCPP
}
// Init
void Init(int MinT, int MaxT, int NSites, int MaxX, double w, double xc, double B0, Rcpp::NumericVector b0, Rcpp::NumericVector b1, Rcpp::NumericVector b2, double pred_a1, double pred_a2, double c, double speed, Rcpp::NumericVector WindAssist, Rcpp::NumericVector WindProb, Rcpp::NumericVector ZStdNorm, Rcpp::NumericVector PStdNorm, Rcpp::NumericVector nTR_x, Rcpp::NumericVector nTR_y, double decError, arma::mat dist, arma::mat bear, double angle, Rcpp::NumericVector y_gain, arma::mat y_expend, Rcpp::NumericVector penalty);
RcppExport SEXP _migrationSDP_Init(SEXP MinTSEXP, SEXP MaxTSEXP, SEXP NSitesSEXP, SEXP MaxXSEXP, SEXP wSEXP, SEXP xcSEXP, SEXP B0SEXP, SEXP b0SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP pred_a1SEXP, SEXP pred_a2SEXP, SEXP cSEXP, SEXP speedSEXP, SEXP WindAssistSEXP, SEXP WindProbSEXP, SEXP ZStdNormSEXP, SEXP PStdNormSEXP, SEXP nTR_xSEXP, SEXP nTR_ySEXP, SEXP decErrorSEXP, SEXP distSEXP, SEXP bearSEXP, SEXP angleSEXP, SEXP y_gainSEXP, SEXP y_expendSEXP, SEXP penaltySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type MinT(MinTSEXP);
    Rcpp::traits::input_parameter< int >::type MaxT(MaxTSEXP);
    Rcpp::traits::input_parameter< int >::type NSites(NSitesSEXP);
    Rcpp::traits::input_parameter< int >::type MaxX(MaxXSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type xc(xcSEXP);
    Rcpp::traits::input_parameter< double >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type pred_a1(pred_a1SEXP);
    Rcpp::traits::input_parameter< double >::type pred_a2(pred_a2SEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type speed(speedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WindAssist(WindAssistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WindProb(WindProbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ZStdNorm(ZStdNormSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type PStdNorm(PStdNormSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nTR_x(nTR_xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nTR_y(nTR_ySEXP);
    Rcpp::traits::input_parameter< double >::type decError(decErrorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bear(bearSEXP);
    Rcpp::traits::input_parameter< double >::type angle(angleSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_gain(y_gainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_expend(y_expendSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type penalty(penaltySEXP);
    Init(MinT, MaxT, NSites, MaxX, w, xc, B0, b0, b1, b2, pred_a1, pred_a2, c, speed, WindAssist, WindProb, ZStdNorm, PStdNorm, nTR_x, nTR_y, decError, dist, bear, angle, y_gain, y_expend, penalty);
    return R_NilValue;
END_RCPP
}
// BackwardIteration
Rcpp::List BackwardIteration();
RcppExport SEXP _migrationSDP_BackwardIteration() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(BackwardIteration());
    return rcpp_result_gen;
END_RCPP
}
// InitSim
void InitSim(int MinT, int MaxT, int NSites, int MaxX, double w, double xc, double B0, Rcpp::NumericVector b0, Rcpp::NumericVector b1, Rcpp::NumericVector b2, double pred_a1, double pred_a2, double c, double speed, Rcpp::NumericVector WindAssist, Rcpp::NumericVector WindProb, Rcpp::NumericVector ZStdNorm, Rcpp::NumericVector PStdNorm, Rcpp::NumericVector nTR_x, Rcpp::NumericVector nTR_y, double decError, arma::mat dist, arma::mat bear, Rcpp::NumericVector y_gain, arma::mat y_expend);
RcppExport SEXP _migrationSDP_InitSim(SEXP MinTSEXP, SEXP MaxTSEXP, SEXP NSitesSEXP, SEXP MaxXSEXP, SEXP wSEXP, SEXP xcSEXP, SEXP B0SEXP, SEXP b0SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP pred_a1SEXP, SEXP pred_a2SEXP, SEXP cSEXP, SEXP speedSEXP, SEXP WindAssistSEXP, SEXP WindProbSEXP, SEXP ZStdNormSEXP, SEXP PStdNormSEXP, SEXP nTR_xSEXP, SEXP nTR_ySEXP, SEXP decErrorSEXP, SEXP distSEXP, SEXP bearSEXP, SEXP y_gainSEXP, SEXP y_expendSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type MinT(MinTSEXP);
    Rcpp::traits::input_parameter< int >::type MaxT(MaxTSEXP);
    Rcpp::traits::input_parameter< int >::type NSites(NSitesSEXP);
    Rcpp::traits::input_parameter< int >::type MaxX(MaxXSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type xc(xcSEXP);
    Rcpp::traits::input_parameter< double >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type pred_a1(pred_a1SEXP);
    Rcpp::traits::input_parameter< double >::type pred_a2(pred_a2SEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type speed(speedSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WindAssist(WindAssistSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type WindProb(WindProbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ZStdNorm(ZStdNormSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type PStdNorm(PStdNormSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nTR_x(nTR_xSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nTR_y(nTR_ySEXP);
    Rcpp::traits::input_parameter< double >::type decError(decErrorSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type bear(bearSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_gain(y_gainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_expend(y_expendSEXP);
    InitSim(MinT, MaxT, NSites, MaxX, w, xc, B0, b0, b1, b2, pred_a1, pred_a2, c, speed, WindAssist, WindProb, ZStdNorm, PStdNorm, nTR_x, nTR_y, decError, dist, bear, y_gain, y_expend);
    return R_NilValue;
END_RCPP
}
// simForaging
arma::vec simForaging(double f_intensity, int time, int site, int x);
RcppExport SEXP _migrationSDP_simForaging(SEXP f_intensitySEXP, SEXP timeSEXP, SEXP siteSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type f_intensity(f_intensitySEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type site(siteSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(simForaging(f_intensity, time, site, x));
    return rcpp_result_gen;
END_RCPP
}
// simFlying
arma::vec simFlying(int decision, int time, int site, int x);
RcppExport SEXP _migrationSDP_simFlying(SEXP decisionSEXP, SEXP timeSEXP, SEXP siteSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type decision(decisionSEXP);
    Rcpp::traits::input_parameter< int >::type time(timeSEXP);
    Rcpp::traits::input_parameter< int >::type site(siteSEXP);
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(simFlying(decision, time, site, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_migrationSDP_Predation", (DL_FUNC) &_migrationSDP_Predation, 5},
    {"_migrationSDP_Init", (DL_FUNC) &_migrationSDP_Init, 27},
    {"_migrationSDP_BackwardIteration", (DL_FUNC) &_migrationSDP_BackwardIteration, 0},
    {"_migrationSDP_InitSim", (DL_FUNC) &_migrationSDP_InitSim, 25},
    {"_migrationSDP_simForaging", (DL_FUNC) &_migrationSDP_simForaging, 4},
    {"_migrationSDP_simFlying", (DL_FUNC) &_migrationSDP_simFlying, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_migrationSDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
