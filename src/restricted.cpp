#include "auxiliary.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <functional>
#include <string>

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(cpp17)]]

// Completely Randomized Design
List crd(NumericVector w, List input) {
  NumericVector probability = w/sum(w);
  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = "CRD",
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Permuted Block Design: PBD (p = b)
List pbd(NumericVector w, List input) {
  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("PBD (")+boost::lexical_cast<std::string>(p)+")";
  double block_size = p*sum(w);

  int j = sum(N) + 1; // current subject's ID
  int k = std::floor((float)((j-1)/block_size));

  NumericVector probability = (w*p*(1+k)-as<NumericVector>(N))/(block_size*(1+k)-(j-1));
  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Block Urn Design: BUD (p = lambda)
List bud(NumericVector w, List input) {
  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("BUD (")+boost::lexical_cast<std::string>(p)+")";

  int j = sum(N) + 1; // current subject's ID
  int k = Rcpp::min(Rcpp::floor(as<NumericVector>(N)/w));

  NumericVector probability = (w*(p+k)-as<NumericVector>(N))/(sum(w)*(p+k)-(j-1));
  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Mass Weight Urn Design: MWUD(p == alpha)
List mwud(NumericVector w, List input) {
  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("MWUD (")+boost::lexical_cast<std::string>(p)+")";

  int j = sum(N) + 1; // current subject's ID

  // target allocation proportions
  NumericVector rho = w/sum(w);

  // tretament groups
  IntegerVector trt = seq_len(rho.size());

  // probabilities of assignment
  NumericVector probability(rho.size());

  for_each(trt.begin(), trt.end(), [p, j, N, rho, &probability](int &k) {
    probability[k-1] = MAX(p*rho[k-1] - N[k-1] + (j-1)*rho[k-1], 0);
  });
  probability = probability/sum(probability);

  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Drop-the-Loser: DL(p == a)
List dl(NumericVector w, List input){
  // processing input
  double p = input["parameter"];    // parameter of the procedure (b)
  IntegerVector N = input["N"];     // current allocation
  IntegerVector urn = input["urn"]; // current urn state

  string procedure = string("DL (")+boost::lexical_cast<std::string>(p)+")";

  // probabilities of assignment
  NumericVector probability(w.size());

  int treatment;

  bool flag = true;
  while(flag){
    treatment = assign_treatment(as<NumericVector>(urn)/sum(urn))-1;

    if (treatment == 0) {
      urn[seq(1, w.size())] = urn[seq(1, w.size())]+(int)(p)*as<IntegerVector>(w);
    }
    else {
      probability = as<NumericVector>(urn)[seq(1, w.size())]/sum(urn[seq(1, w.size())]);
      urn[treatment] = urn[treatment]-1;
      flag = false;
    }
  }

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment,
    _["urn"] = urn
  );
}



// Generalized Drop-the-Loser: GDL(p == a)
List gdl(NumericVector w, List input){
  // processing input
  double p = input["parameter"];    // parameter of the procedure (b)
  IntegerVector N = input["N"];     // current allocation
  IntegerVector urn = input["urn"]; // current urn state

  string procedure = string("GDL (")+boost::lexical_cast<std::string>(p)+")";

  // probabilities of assignment
  NumericVector probability(w.size());

  int treatment;

  bool flag = true;
  while(flag){
    // urn update (if there are urns with negative number of balls)
    for(int u = 0; u < urn.size();u++) { urn[u] = MAX(urn[u], 0); }

    treatment = assign_treatment(as<NumericVector>(urn)/sum(urn))-1;

    if (treatment == 0) {
      urn[seq(1, w.size())] = urn[seq(1, w.size())]+(int)(p)*as<IntegerVector>(w);
    }
    else {
      probability = as<NumericVector>(urn)[seq(1, w.size())]/sum(urn[seq(1, w.size())]);
      urn[treatment] = urn[treatment]-1;
      flag = false;
    }
  }

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment,
    _["urn"] = urn
  );
}


// Doubly-Adaptive Biased Coind Design: DBCD(p == gamma)
List dbcd(NumericVector w, List input) {
  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("DBCD (")+boost::lexical_cast<std::string>(p)+")";

  int j = sum(N) + 1; // current subject's ID

  // target allocation proportions
  NumericVector rho = w/sum(w);

  // tretament groups
  IntegerVector trt = seq_len(rho.size());

  // probabilities of assignment
  NumericVector probability(rho.size());

  if (is_true(all(N > 0))) {
    probability = rho*Rcpp::pow((rho/(as<NumericVector>(N)/(j-1))), p);
    probability = probability/sum(probability);
  }
  else {
    probability = rho;
  }

  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Minimum Quadratic Distance: MinQD(p == eta)
List min_qd(NumericVector w, List input){
  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  if (p < 0 || p > 1) {
    throw invalid_argument(
        string("The value of randomizaton procedure parameter must belong to [0, 1]: ") +
          " current value equals to " + boost::lexical_cast<std::string>(p)
    );
  }
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("MinQD (")+boost::lexical_cast<std::string>(p)+")";

  int j = sum(N) + 1; // current subject's ID

  // target allocation proportions
  NumericVector rho = w/sum(w);

  // the hypothetical "lack of balance"
  NumericVector B(rho.size());

  // tretament groups
  IntegerVector trt = seq_len(rho.size());

  // probabilities of assignment
  NumericVector probability(rho.size());

  for_each(trt.begin(), trt.end(), [j, N, rho, &B](int &k){
    IntegerVector N1 = N;
    N1[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
    N1[k-1] -= 1;
  });

  double mu;

  if (var(B) <= 1e-16) {
    probability = rho;
  }
  else {
    mu = 2/((w.size()-1)*var(B))*p*(sum(B*rho)-min(B));
    probability = rho-0.5*mu*(B-mean(B));
  }
  if (!is_true(all(probability >= 0))) {

    LogicalVector id = (probability >= 0);
    probability[!id] = 0;
    double Q = sum(as<NumericVector>(rho[id])) - 1;
    double Kw = sum(id);
    double Bw = sum(as<NumericVector>(B[id]));

    mu = 2*(p*(min(B)-sum(B*rho)) + Q/Kw*sum(B))/sum(B*(B-Bw/Kw));
    for_each(trt.begin(), trt.end(), [id, rho, B, mu, Q, Kw, Bw, &probability](int &k){
      if (id[k-1]) {
        probability[k-1] = rho[k-1]-Q/Kw+mu/2*(B[k-1]-Bw/Kw);
      }
    });

  }

  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


// Maximum Entropy: MaxEnt(p == eta)
List max_ent(NumericVector w, List input){

  // processing input
  double p = input["parameter"]; // parameter of the procedure (b)
  if (p < 0 || p > 1) {
    throw invalid_argument(
        string("The value of randomizaton procedure parameter must belong to [0, 1]: ") +
          " current value equals to " + boost::lexical_cast<std::string>(p)
    );
  }
  IntegerVector N = input["N"];  // current allocation

  string procedure = string("MaxEnt (")+boost::lexical_cast<std::string>(p)+")";

  int j = sum(N) + 1; // current subject's ID

  // target allocation proportions
  NumericVector rho = w/sum(w);

  // the hypothetical "lack of balance"
  NumericVector B(rho.size());

  // tretament groups
  IntegerVector trt = seq_len(rho.size());

  // probabilities of assignment
  NumericVector probability(rho.size());

  for_each(trt.begin(), trt.end(), [j, N, rho, &B](int &k){
    IntegerVector N1 = N;
    N1[k-1] += 1;
    B[k-1] = max(Rcpp::abs(as<NumericVector>(N1)/(double)(j)-rho));
    N1[k-1] -= 1;
  });

  if (p != 1) {
    if (var(B) <= 1e-16) {
      probability = rho;
    }
    else {
      // we have to find a zero of a one-variable (mu) function.
      // bisection method is used -- implemeted in the file algorithms.cpp

      // function to find zero of
      auto fcn = [p, B, rho](double mu) {
        return min(B)*p + (1-p)*sum(rho*B) - sum(B*rho*pow(exp(B), -mu))/sum(rho*pow(exp(B), -mu));
      };

      auto der = [p, B, rho](double mu) {
        return sum(B*B*rho*pow(exp(B), -mu))/sum(rho*pow(exp(B), -mu)) -
          pow(sum(B*rho*pow(exp(B), -mu))/sum(rho*pow(exp(B), -mu)), 2);
      };

      double mu = newton_raphson(10, fcn, der, 1e-5, 1e-5, 50);
      probability = rho*pow(exp(B), -mu)/sum(rho*pow(exp(B), -mu));
    }
  }
  else {
    probability[B == min(B)] = rho[B == min(B)];
    probability = probability/sum(probability);
  }

  int treatment = assign_treatment(probability);

  return List::create(
    _["procedure"] = procedure,
    _["probability"] = probability,
    _["treatment"] = treatment
  );
}


//[[Rcpp::export(name=.restricted_one_simulation_cpp)]]
List restricted_one_simulation(
    NumericVector w,
    int nsbj,
    std::string procedure,
    double parameter) {

  CharacterVector rnd_procedure = CharacterVector::create(
    "CRD",
    "PBD",
    "BUD",
    "MWUD",
    "DL",
    "GDL",
    "DBCD",
    "MinQD",
    "MaxEnt"
  );


  // define a randomization function
  typedef boost::function< List(NumericVector, List) > RandomizationFcn;
  std::vector< RandomizationFcn > rnd_functions;
  rnd_functions.push_back(boost::bind(&crd, _1, _2));
  rnd_functions.push_back(boost::bind(&pbd, _1, _2));
  rnd_functions.push_back(boost::bind(&bud, _1, _2));
  rnd_functions.push_back(boost::bind(&mwud, _1, _2));
  rnd_functions.push_back(boost::bind(&dl, _1, _2));
  rnd_functions.push_back(boost::bind(&gdl, _1, _2));
  rnd_functions.push_back(boost::bind(&dbcd, _1, _2));
  rnd_functions.push_back(boost::bind(&min_qd, _1, _2));
  rnd_functions.push_back(boost::bind(&max_ent, _1, _2));

  RandomizationFcn randomization_fcn;

  for(int rnd = 1; rnd <= rnd_procedure.size(); rnd++) {
    if (procedure == as<std::string>(rnd_procedure[rnd-1])) {
      randomization_fcn = rnd_functions[rnd-1];
      break;
    }
    if (rnd == 9) {
      throw invalid_argument(
          string("The value of input parameter procedure must be one of ") +
            "c('CRD', 'PBD', 'BUD', 'MWUD', 'DL', 'GDL', 'DBCD', 'MinQD', 'MaxEnt'): " +
            "current value is '" + procedure + "'"
      );
    }
  }


  // number of treatments
  int ntrt = w.size();

  // target allocation proportions
  NumericVector rho = w/sum(w);

  // allocations at each steps
  IntegerMatrix N(nsbj+1, ntrt);

  // probabilities of assignment at each step
  NumericMatrix probability(nsbj, ntrt);

  // assigned treatment at each step
  IntegerVector treatment(nsbj);

  // forcing index at each step
  NumericVector fi_(nsbj);
  NumericVector fi(nsbj);

  // imbalance
  NumericVector imbalance(nsbj);

  // initial urn state (used by urn-based procedues, e.g. DL, GDL)
  NumericVector urn1(ntrt+1);
  urn1[0] = 1;
  urn1[seq(1, w.size())] = w;

  // run simulation
  IntegerVector subject = seq_len(nsbj);

  List input, out;
  for (int j = 1; j <= subject.size(); j++) {
    IntegerVector Nj = N(j-1,_); // a vector of current allocation
    input = List::create(
      _["parameter"] = parameter,
      _["N"] = Nj,
      _["urn"] = urn1
    );

    out = randomization_fcn(w, input);
    treatment[j-1] = out["treatment"];
    probability(j-1, _) = as<NumericVector>(out["probability"]);
    Nj[treatment[j-1]-1] += 1;
    N(j,_) = Nj;

    // operational characteristics
    fi_[j-1] = sum(Rcpp::pow(probability(j-1,_)-rho, 2));
    fi[j-1] = mean(fi_[seq(0, j-1)]);

    imbalance[j-1] = sqrt((float)sum(Rcpp::pow(as<NumericVector>(Nj) - j*rho, 2)));

    if (procedure == "DL" || procedure == "GDL"){
      urn1 = out["urn"];
    }
  }

  List output = List::create(
    _["procedure"] = out["procedure"],
    _["subject"] = subject,
    _["treatment"] = treatment,
    _["FI"] = fi,
    _["imbalance"] = imbalance,
    _["allocation"] = N(seq(1, nsbj), _),
    _["probability"] = probability
  );

  return output;
}


//[[Rcpp::export(name=.restricted_multiple_simulations_cpp)]]
List restricted_multiple_simulations(
    NumericVector w,
    int nsbj,
    std::string procedure,
    double parameter,
    int nsim) {

  List out(nsim);

  IntegerVector sim = seq_len(nsim);
  for_each(sim.begin(), sim.end(), [&out, w, nsbj, procedure, parameter](int &s){
    out[s-1] = restricted_one_simulation(w, nsbj, procedure, parameter);
  });

  return out;
}
