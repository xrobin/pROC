#include <Rcpp.h> 
using namespace Rcpp;

double MWkernel(double x, double y) {
  if (y == x) {return(0.5);}
	else if (y < x) {return(1);}
	else {return(0);}
}


// [[Rcpp::export]]
List delongPlacementsCpp(List roc) {  
  NumericVector cases = roc["cases"];
	NumericVector controls = roc["controls"];
	std::size_t xsize = cases.size();
	std::size_t ysize = controls.size();
	double sum = 0;
	double MW;
	std::vector<double> X(xsize, 0.0);
	std::vector<double> Y(ysize, 0.0);

	for (std::size_t i = 0; i < xsize; i++) {
	  for (std::size_t j = 0; j < ysize; j++) {
      MW = MWkernel(cases[i], controls[j]);
      X[i] += MW ;
	    Y[j] += MW;
	    sum += MW;
	  }
	}
  for (std::size_t i = 0; i < xsize; i++) {
    X[i] /= ysize;
  }
	for (std::size_t j = 0; j < ysize; j++) {
    Y[j] /= xsize;
	}
	List ret;
	ret["theta"] = sum / (xsize * ysize);
	ret["X"] = X;
	ret["Y"] = Y;
	return(ret);
}
