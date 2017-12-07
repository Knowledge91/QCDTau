// Copyright 2017 Dirk Hornung

#include "./chisquared.h"

typedef Weights W;

using ROOT::Math::Minimizer;
using ROOT::Math::Factory;
using ROOT::Math::Functor;
using std::cout;
using std::endl;
using std::vector;
using std::function;

// chisquared test
double RosenBrock(const double *xx ) {
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

Chisquared::Chisquared(Constants constants, vector<double> s0Set) :
    constants_(constants), s0_set_(s0Set), jacobian_matrix_(82, 9),
    covariance_matrix_(9, 9), data_(ExperimentalData()),
    experiment_(ExperimentalMoment(3.1570893124000001, W::WD00)),
    inverse_covariance_matrix_(9, 9) {
  FillJacobianMatrix();
  FillCovarianceMatrix();

  // remove correlations with R_tau, V+A in Aleph fit
  for (int i = 0; i < 9; i++) {
    covariance_matrix_(1, i) = 0.;
    covariance_matrix_(i, 1) = 0.;
  }

  cout << covariance_matrix_ << endl;
  Numerics::InvertMatrix(covariance_matrix_, inverse_covariance_matrix_);
  cout << inverse_covariance_matrix_ << endl;
}

void Chisquared::Minuit() {
  Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");

  // set tolerances
  min->SetMaxFunctionCalls(1000000);  // for Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // function wrapper
  Functor f(&RosenBrock, 2);
  double step[2] = {0.01, 0.01};
  // starting point
  double variable[2] = {-1., 1.2};

  min->SetFunction(f);

  // set free variables to be minimized
  min->SetVariable(0, "x", variable[0], step[0]);
  min->SetVariable(1, "y", variable[1], step[1]);

  // minimize!
  min->Minimize();

  const double *xs = min->X();
  cout << "Minimum f(" << xs[0] << "," << xs[1] << "): "
            << min->MinValue() << endl;
}

void Chisquared::PrintS0Set() {
    cout << "S0Set: " << GetS0SetSize() << " points" << endl;
    for (int i = 0; i < GetS0SetSize(); i++) {
      cout << "S0(" << i+1 << ") \t" << s0_set_[i] << endl;
    }
  }

void Chisquared::FillJacobianMatrix() {
  function<double(double)> weight = W::WD00;
  for (int i = 0; i < GetS0SetSize(); i++) {
    double s0 = s0_set_[i];
    ExperimentalMoment experimentalMoment(s0, weight);
    vector<double> jacobian_vector =
        experimentalMoment.GetJacobianVector();
    for (int j = 0; j < data_.GetNumberOfDataPoints()+2; j++) {
      jacobian_matrix_(j, i) = jacobian_vector[j];
    }
  }
}

void Chisquared::FillCovarianceMatrix() {
  for (int i = 0; i < GetS0SetSize(); i++) {
    for (int j = 0; j < GetS0SetSize(); j++) {
      for (int k = 0; k < data_.GetNumberOfDataPoints()+2; k++) {
        for (int l = 0; l < data_.GetNumberOfDataPoints()+2; l++) {
          covariance_matrix_(i, j) += jacobian_matrix_(k, i)*
              data_.GetErrorMatrix(k, l)*jacobian_matrix_(l, j);
        }
      }
    }
  }
}
