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


namespace Ublas {

  using boost::numeric::ublas::matrix;

class Chisquared {
 public:
  Chisquared(Constants constants) : constants_(constants),
      experiment_(ExperimentalMoment()) {}

 private:
  double invertedCovarianceMatrix[80][80];
  Constants constants_;
  ExperimentalMoment experiment_;
};

}  // namespace Ublas

#endif
