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
#pragma once 

#include <algorithm> // std::iota, std::sort
#include <string>
#include <vector>

namespace pROC {
	template <class P> std::vector<int> getPredictorOrder(const P& predictor, const std::string& direction) {
		std::vector<int> index(predictor.nTotal);
		std::iota(index.begin(), index.end(), 0);
		if (direction == ">") {
			std::sort(index.begin(), index.end(), PredictorComparator<P>(predictor));
		}
		else {
			std::sort(index.begin(), index.end(), PredictorReverseComparator<P>(predictor));
		}
		
		return index;
	}
}
