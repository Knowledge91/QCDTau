// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
#define PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_

#include <vector>
#include "Constants.h"
#include "Weights.h"
#include "./Numerics.h"
#include "./experimental_data.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>



typedef Constants C;
typedef Weights W;
typedef Numerics N;

namespace Ublas {

  using boost::numeric::ublas::matrix;
  using std::vector;
  using std::cout;

class ExperimentalMoment {
 public:
  ExperimentalMoment(double s0, function<double(double)> weight) : s0_(s0),
      weight_(weight), data_(ExperimentalData()) {}

  // get Spectral-moment for -s0, -weight(x)
  double SpectralMoment(const double s0, function<double(double)> weight) {
    int N = getBinNumber(s0);
    double sum = 0;
    for (int i=0; i <= N; i++) {
      sum += wRatio(s0, weight, i)*data_.GetSfm2(i);
    }
    return C::sTau/s0/C::Be*sum;
  }

  template <typename Func>
  double wRatio(double s0, Func weight, int i) {
    return
        s0/ C::sTau*
        (weight((data_.GetSbin(i)-data_.GetDSbin(i)/2.)/s0)
         -weight((data_.GetSbin(i)+data_.GetDSbin(i)/2.)/s0) )
        /
        (W::WD00((data_.GetSbin(i)-data_.GetDSbin(i)/2.)/C::sTau)
         - W::WD00((data_.GetSbin(i)+data_.GetDSbin(i)/2.)/C::sTau));
  }

  double s0_;
  function<double(double)> weight_;
  ExperimentalData data_;

  // get last included bin from s0
  int getBinNumber(double s0) {
    double pos = 0;
    for (int i=0; i < 80; i++) {
      pos += data_.GetDSbin(i);
      if (pos >= s0) {
        return i+1;
      }
    }
    return 80;
  }

  vector<double> GetJacobianVector() {
    vector<double> jacobian(80);
    for (int i = 0; i < 80; i++) {
      int N = getBinNumber(s0_);
      if (i < N) {
        jacobian[i] = Constants::sTau/Constants::Be/s0_*wRatio(s0_, weight_, i);
      } else {
        jacobian[i] = 0.;
      }
    }
    return jacobian;
  }
};

}  // namespace Ublas

using Ublas::ExperimentalMoment;

#endif  // PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
