#pragma once

#include <Rcpp.h>
using namespace Rcpp;

#include "auc.h"

namespace Rcpp {
    template <> AucParams as(SEXP someParams);
}
