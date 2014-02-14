#include <vector>
#include <utility>

/* class BootstrapGenerator {
#include <random>
  std::mt19937 rng = std::mt19937(std::random_device{}());
  std::uniform_int_distribution<size_t> uniform_dist;
  std::vector<double> myVector;
  
  public:
    //BootstrapGenerator() : rd(time(0)), rng(rd()) {};
  BootstrapGenerator(const std::vector<double>& aVector) : uniform_dist(0, aVector.size() - 1), myVector(aVector) {}
  double operator()() { return myVector[uniform_dist(rng)]; }
};*/



void setRandomSample(const std::vector<double>& from, std::vector<double>& to);
void setRandomUnpairedSample(const std::vector<double>& fromControls, const std::vector<double>& fromCases,
                             std::vector<double>& toControls, std::vector<double>& toCases);

std::vector<double> computeThresholds(const std::vector<double>& controls, const std::vector<double>& cases);

std::pair<std::vector<double>, std::vector<double>> computeSeSp(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

double computeAuc(std::pair<std::vector<double>, std::vector<double>>);