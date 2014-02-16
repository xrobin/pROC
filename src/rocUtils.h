#include <vector>
#include <utility>


std::vector<double> computeThresholds(const std::vector<double>& predictor);
std::vector<double> computeThresholds(const std::vector<double>& controls, const std::vector<double>& cases);
std::vector<double> computeThresholds(const std::vector<double>& controls, const std::vector<double>& cases, 
                                      const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx);

std::pair<std::vector<double>, std::vector<double>> computeSeSp(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
std::pair<std::vector<double>, std::vector<double>> computeSeSp(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, 
                                                                const std::vector<size_t>& controlsIdx, const std::vector<size_t>& casesIdx);
