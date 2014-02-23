/* pROC: Tools Receiver operating characteristic (ROC curves) with
   (partial) area under the curve, confidence intervals and comparison. 
   Copyright (C) 2014 Xavier Robin.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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