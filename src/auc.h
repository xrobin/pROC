#include <Rcpp.h>
#include <vector>

struct AucParams {
  bool partial, focusOnSp, correct;
  double from, to;
  explicit AucParams(const Rcpp::List& l);
};


double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, const AucParams& aucParams);
double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, 
             const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx, 
             const AucParams& aucParams);

double computeAuc(const std::pair<std::vector<double>, std::vector<double>>&, const AucParams& aucParams);