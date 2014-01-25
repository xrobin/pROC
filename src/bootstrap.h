#pragma once

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector bootstrapAucStratified(Rcpp::NumericVector controls, Rcpp::NumericVector cases, int bootN);