#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List rocUtilsPerfsAllC(NumericVector thresholds, NumericVector controls, NumericVector cases, std::string direction) {
  NumericVector se(thresholds.size());
  NumericVector sp(thresholds.size());
  long tp, tn;
  long i; // iterator over cases & controls
  
  if (direction == ">") {
    for (long t = 0; t < thresholds.size(); t++) {
      double threshold = thresholds(t);
        tp = 0;
        for (i = 0; i < cases.size(); i++) {
          if (cases(i) <= threshold) {
            tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) > threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
 else {
    for (long t = 0; t < thresholds.size(); t++) {
      double threshold = thresholds(t);
      tp = 0;
      for (i = 0; i < cases.size(); i++) {
        if (cases(i) >= threshold) {
          tp++;
        }
      }
      se(t) = (double)tp / cases.size();
      long tn = 0;
      for (i = 0; i < controls.size(); i++) {
        if (controls(i) < threshold) {
          tn++;
        }
      }
      sp(t) = (double)tn / controls.size();
    }
  }
  List ret;
  ret["se"] = se;
  ret["sp"] = sp;
  return(ret);
}