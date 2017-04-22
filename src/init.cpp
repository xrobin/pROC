/* pROC: Tools Receiver operating characteristic (ROC curves) with
 (partial) area under the curve, confidence intervals and comparison. 
 Copyright (C) 2014 Xavier Robin
 
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

#include <R.h>
#include <R_ext/Rdynload.h>

/** This function is required in R >= 3.4 to avoid a NOTE in 
 * R CMD check.
 * Source:
 * https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661
 */
void R_init_pROC(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}