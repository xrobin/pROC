#include <Rcpp.h>
using namespace Rcpp ;

#include "Predictor.h"

// [[Rcpp::export]]
NumericVector runit_Predictor_bracketOperator(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const int i) {
	Predictor aPredictor(someControls, someCases);
	return aPredictor[i];
}

// [[Rcpp::export]]
NumericVector runit_Predictor_bracketOperatorVector(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases) {
	Predictor aPredictor(someControls, someCases);
	std::vector<double> ret;
	size_t maxI = someControls.size() + someCases.size();
	ret.reserve(maxI);
	for (size_t i = 0; i < maxI; ++i) {
		ret.push_back(aPredictor[i]);
	}
	return wrap(ret);
}

// [[Rcpp::export]]
IntegerVector runit_Predictor_getOrder(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases) {
	Predictor aPredictor(someControls, someCases);
	return wrap(aPredictor.getOrder());
}

// [[Rcpp::export]]
NumericVector runit_Predictor_getControls(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases) {
	Predictor aPredictor(someControls, someCases);
	return wrap(aPredictor.getControls());
}

// [[Rcpp::export]]
NumericVector runit_Predictor_getCases(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases) {
	Predictor aPredictor(someControls, someCases);
	return wrap(aPredictor.getCases());
}

// [[Rcpp::export]]
bool runit_Predictor_isControl(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const int i) {
	Predictor aPredictor(someControls, someCases);
	return aPredictor.isControl(i);
}

// [[Rcpp::export]]
bool runit_Predictor_isCase(const Rcpp::NumericVector& someControls, const Rcpp::NumericVector& someCases, const int i) {
	Predictor aPredictor(someControls, someCases);
	return aPredictor.isCase(i);
}