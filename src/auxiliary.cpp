#include "auxiliary.hpp"

// function to assign a treatment, given probabilities
int assign_treatment(NumericVector probability) {
  int k; // treatment id
  NumericVector cumulative_prob(probability.size()+1);

  double u = runif(1)[0];
  for(k = 1; k <= probability.size(); k++) {
    cumulative_prob[k] = cumulative_prob[k-1] + probability[k-1];
    if (cumulative_prob[k-1] < u && u <= cumulative_prob[k]) {
      break;
    }
  }
  return k;
}


// function to perform ANOVA test
//[[Rcpp::export(name=.anova_test_cpp)]]
int anova_test(
    IntegerVector treatment,
    NumericVector response,
    int K,
    double alpha) {

  // processing input
  int N = treatment.size(); // sample size

  IntegerVector n(K);
  NumericVector Y_(K);
  double Y = mean(response);

  double MSw = 0;
  double MSb = 0;

  int reject = R_NaInt;

  if (N > K) {
    for(int k = 1; k <= K; k++) {
      NumericVector Yk = response[treatment == k];
      n[k-1] = Yk.size();
      if (n[k-1] <= 1) break;
      else {
        Y_[k-1] = mean(Yk);
        MSw += sum(pow(Yk - Y_[k-1], 2));
        MSb += n[k-1]*pow(Y-Y_[k-1], 2);
      }
    }
    double F = (MSb/(K-1))/(MSw/(N-K));
    reject = (int)(F > Rcpp::qf(NumericVector::create(1-alpha), K-1, N-K))[0];
  }

  return reject;
}



// function to add responses, given treatments, and to pervorm ANOVA test

//[[Rcpp::export(name=.anova_test_per_subject_cpp)]]
List anova_test_per_subject(
    List treatment_list,
    List response_list,
    int K,
    double alpha,
    bool time_drift
) {

  IntegerMatrix reject_mat(as<IntegerVector>(treatment_list[0]).size(), treatment_list.size());
  IntegerVector sim = seq_len(treatment_list.size());
  for_each(sim.begin(), sim.end(), [&reject_mat, treatment_list, response_list, K, alpha, time_drift](int &s) {
    IntegerVector treatment = treatment_list[s-1];
    NumericMatrix response_mat = response_list[s-1];
    NumericVector response(treatment.size());

    // perform ANOVA test at each allocation step
    IntegerVector subject = seq_len(treatment.size());
    IntegerVector reject(subject.size());
    for_each(subject.begin(), subject.end(), [treatment, response_mat, &response](int &sbj){
      response[sbj-1] = response_mat(sbj-1, treatment[sbj-1]-1);
    });
    if (time_drift) {
      response += as<NumericVector>(subject)/subject.size();
    }
    for_each(subject.begin(), subject.end(), [&reject, treatment, response, K, alpha](int &sbj){
      reject[sbj-1] = anova_test(treatment[seq(0, sbj-1)], response[seq(0, sbj-1)], K, alpha);
    });
    reject_mat.column(s-1) = reject;
  });

  IntegerVector subject = seq_len(reject_mat.nrow());
  NumericVector error(subject.size());
  for_each(subject.begin(), subject.end(), [&error, reject_mat](int &sbj){
    error[sbj-1] = mean(reject_mat.row(sbj-1));
  });

  return List::create(
    _["subject"] = subject,
    _["error"] = error
  );
}

// function to solve nonlinear equation
double newton_raphson(
    double x0,
    std::function<double (double)> fcn,
    std::function<double (double)> der,
    double xtol,
    double ftol,
    int maxiter){

  double a = 0.05;
  double x = x0-a*fcn(x0)/der(x0);
  int iter = 0;
  while ((abs(fcn(x)) > ftol) && (abs(x-x0) > xtol) && (iter <= maxiter)) {
    x0 = x;
    x = x0-a*fcn(x0)/der(x0);
    iter += 1;
  }

  return x;
}
