#include <Rcpp.h>
#include <vector>


double aucCC(const std::vector<double>& controls, const std::vector<double>& cases);
double aucCC(const std::vector<double>& controls, const std::vector<double>& cases, 
             const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx);

double computeAuc(const std::pair<std::vector<double>, std::vector<double>>&);