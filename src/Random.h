#include <vector>


void setRandomIdx(const size_t size, std::vector<size_t>& idxVector);
void setRandomNonStratifiedSample(const size_t controlsSize, const size_t casesSize,
                             std::vector<size_t>& controlsIdx, std::vector<size_t>& casesIdx);