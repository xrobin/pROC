#include <vector>
#include <string>
#include "Predictor.h"

using std::vector;
using std::string;

vector<size_t> Predictor::getOrder() const {
    vector<size_t> index(nTotal);
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), PredictorComparator(*this));
    return index;
}

vector<size_t> Predictor::getOrder(const string direction) const {
    vector<size_t> index(nTotal);
    std::iota(index.begin(), index.end(), 0);
    if (direction == ">") {
      std::sort(index.begin(), index.end(), PredictorComparator(*this));
    }
    else {
      std::sort(index.begin(), index.end(), PredictorReverseComparator(*this));
    }
    
    return index;
}