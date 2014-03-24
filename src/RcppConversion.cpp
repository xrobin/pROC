#include <Rcpp.h>
using namespace Rcpp;

#include <string>
using std::string;

#include "RcppConversions.h"
#include "auc.h"

namespace Rcpp {
    template <> AucParams as(SEXP rParams) {
    	AucParams cppParams;
    	List listParams = as<List>(rParams);

		try {
			bool test = listParams["partial.auc"];
			if (! test) {
				cppParams.partial = false;
				return cppParams;
			}
			else {
				stop("Reached a line that should be unreachable (partial.auc=TRUE). Please report this bug to the maintainer of pROC. Type packageDescription(\"pROC\", fields=\"Maintainer\") to obtain this information.");
			}
		}
		catch (const Rcpp::not_compatible&) { // the conversion to bool will throw for partial auc where test is a numeric vector of length 2.
			cppParams.partial = true;
			NumericVector test = listParams["partial.auc"];
			cppParams.from = test[0];
			cppParams.to = test[1];
			
			cppParams.correct = listParams["partial.auc.correct"];
			string paucFocus = listParams["partial.auc.focus"];
			if (paucFocus == "specificity") {
				cppParams.focusOnSp = true;
			}
			else if (paucFocus == "sensitivity") {
				cppParams.focusOnSp = false;
			}
			else {
				string errMsg = string("Invalid partial.auc.focus: ") + paucFocus + ". This probably denotes a bug in pROC. If so, please type packageDescription(\"pROC\", fields=\"Maintainer\") to and report it to the maintainer.";
				stop(errMsg);
			}
			return cppParams;
		}
		stop("Reached a line that should be unreachable (end of as<AucParams>). Please report this bug to the maintainer of pROC. Type packageDescription(\"pROC\", fields=\"Maintainer\") to obtain this information.");
    }
}   
