// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_CHISQUARED_H_
#define PROGRAM_SRC_CHISQUARED_H_

#include <vector>
#include "Numerics.h"
#include "./experimental_moment.h"
#include "./theoretical_moment.h"

#include "storage_adaptors.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

double RosenBrock(const double *xx ) {
  const Double_t x = xx[0];
  const Double_t y = xx[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}


namespace ROOT::Math {

  using boost::numeric::ublas::matrix;

class Chisquared {
 public:
  explicit Chisquared(Constants constants) : constants_(constants),
    experiment_(ExperimentalMoment()) {}

  static void minuit() {
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
    std::cout << "Minimum f(" << xs[0] << "," << xs[1] << "): "
              << min->MinValue() << std::endl;
  }

 private:
  double invertedCovarianceMatrix[80][80];
  Constants constants_;
  ExperimentalMoment experiment_;
};


}  // namespace ROOT::Math

using ROOT::Math::Chisquared;

#endif  // PROGRAM_SRC_CHISQUARED_H_
