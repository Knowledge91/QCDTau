// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
#define PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_

#include <vector>
#include <functional>
#include "./constants.h"
#include "./weights.h"
#include "./experimental_data.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace experimental_moment_wrapper {
  typedef Constants C;
  typedef Weights W;

  using std::function;
  using std::vector;

class ExperimentalMoment {
 public:
  ExperimentalMoment(double s0, function<double(double)> weight);

  // get Spectral-moment for s0_ and weight_
  double GetSpectralMoment();

  // calculates the Jacobian vector for the Gauss-Error-Propagation
  vector<double> GetJacobianVector();

 private:
  // calculates the needed weight_ ratio for the spectral moment
  //
  // in general we use as numerator W::WD00 ( weight_ / WD00 )
  double wRatio(int i);

  // get last included bin from s0_
  //
  // The experimental data is discrete. As we are summing over the data until
  // s0_ we need to know when to stop summing.
  int getBinNumber();

  double s0_;
  std::function<double(double)> weight_;
  ExperimentalData data_;
};

}  // namespace experimental_moment_wrapper
using experimental_moment_wrapper::ExperimentalMoment;

#endif  // PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
