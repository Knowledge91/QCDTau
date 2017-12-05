// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_CHISQUARED_H_
#define PROGRAM_SRC_CHISQUARED_H_

#include <vector>
#include "./numerics.h"
#include "./weights.h"
#include "./experimental_moment.h"
#include "./theoretical_moment.h"
#include "./experimental_data.h"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"




class Chisquared {
 public:
  Chisquared(Constants constants, std::vector<double> s0Set);

  double GetS0SetSize() const { return s0_set_.size(); }

  double GetJacobianMatrix(int i, int j) const {
    return jacobian_matrix_(i, j);
  }

  double GetCovarianceMatrix(int i, int j) const {
    return covariance_matrix_(i, j);
  }

  double GetInvertedCovarianceMatrix(int i, int j) const {
    return inverted_covariance_matrix_(i, j);
  }

  void PrintS0Set();

  static void Minuit();

  ExperimentalData data_;
  std::vector<double> s0_set_;
  boost::numeric::ublas::matrix<double> jacobian_matrix_;
  boost::numeric::ublas::matrix<double> covariance_matrix_;
  boost::numeric::ublas::matrix<double> inverted_covariance_matrix_;
  Constants constants_;
  ExperimentalMoment experiment_;

  void FillJacobianMatrix();

  void FillCovarianceMatrix();
};

#endif  // PROGRAM_SRC_CHISQUARED_H_
