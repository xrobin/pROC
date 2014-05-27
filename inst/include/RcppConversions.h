#pragma once

#include <Rcpp.h>

#include <pROC/auc.h>
#include <pROC/ROC.h>
#include <pROC/Predictor.h>

namespace Rcpp {
	template <> pROC::AucParams as(SEXP someParams);

	template <> pROC::ROC<> as(SEXP aROC);
	template <> SEXP wrap(const pROC::ROC<> &aROC);
}
