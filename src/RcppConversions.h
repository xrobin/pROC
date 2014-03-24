#pragma once

#include <Rcpp.h>
using namespace Rcpp;

#include "auc.h"
#include "ROC.h"
#include "Predictor.h"

namespace Rcpp {
	template <> AucParams as(SEXP someParams);

	template <> ROC<Predictor> as(SEXP aROC);
	template <> SEXP wrap(const ROC<Predictor> &aROC);
}
