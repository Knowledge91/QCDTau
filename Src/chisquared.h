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



class Chisquared {
 public:
  Chisquared(Constants constants) : constants_(constants), experiment_(ExperimentalMoment()) {};

  void build(double s0) {
    MatDoub_IO covarianceMatrix = experiment_.covarianceMatrix;
    MatDoub_IO invertedCovarianceMatrix = invertCovarianceMatrix(covarianceMatrix);
  }
  
  
 private:
  double invertedCovarianceMatrix[80][80];
  Constants constants_;
  ExperimentalMoment experiment_;

  MatDoub_IO invertCovarianceMatrix(MatDoub_IO &covarianceMatrix) {
    MatDoub_IO invertedCovarianceMatrix;
    Numerics::gaussj(covarianceMatrix, invertedCovarianceMatrix);
    cout << "covarianceMatrix[0][0] \t" << covarianceMatrix[0][0] << endl;
    cout << "invertedCovarianceMatrix[0][0] \t" << invertedCovarianceMatrix[0][0] << endl;
    return invertedCovarianceMatrix;
  }
};

#endif
