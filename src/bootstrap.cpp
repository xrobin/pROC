#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <random>

using std::vector;
using namespace Rcpp;

class BootstrapGenerator {
  std::mt19937 rng = std::mt19937(std::random_device{}());
  std::uniform_int_distribution<size_t> uniform_dist;
  std::vector<double> vector;
  
  public:
    //BootstrapGenerator() : rd(time(0)), rng(rd()) {};
    BootstrapGenerator(const std::vector<double>& vector) : uniform_dist(0, vector.size() - 1), vector(vector) {}
    double operator()() { return vector[uniform_dist(rng)]; }
};


// [[Rcpp::export]]
std::vector<double> bootstrapAucStratified(const std::vector<double>& controls, const std::vector<double>& cases, const int bootN) {
  // keep all AUCs in a vector of size bootN
  vector<double> aucs;
  aucs.reserve(bootN);
  
  std::vector<double> sampleControls(controls.size()), 
                      sampleCases(cases.size());
  
  // Generate sampling vector
  BootstrapGenerator controlsGenerator(controls),
                     casesGenerator(cases);

  for (int i = 0; i < bootN; i++) {
    std::generate(sampleControls.begin(), sampleControls.end(), controlsGenerator);
    std::generate(sampleCases.begin(), sampleCases.end(), casesGenerator);
  
    // Compute SE/SP of sample
    
    // Compute AUC
  }
  
  return aucs;
}