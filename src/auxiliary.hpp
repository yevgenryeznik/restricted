#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// function to assign a treatment, given probabilities
int assign_treatment(NumericVector probability);


// function to perform ANOVA test
int anova_test(
    IntegerVector,
    NumericVector,
    int,
    double);


// function to add responses, given treatments, and to pervorm ANOVA test
List anova_test_per_subject(
  List,
  List,
  int,
  double,
  bool);


// function to solve nonlinear equation
double newton_raphson(
    double,
    std::function<double (double)>,
    std::function<double (double)>,
    double,
    double,
    int);
