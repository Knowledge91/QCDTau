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


typedef Constants C;
typedef Weights W;


class ExperimentalMoment {
 public:
  ExperimentalMoment(double s0, std::function<double(double)> weight);

  // get Spectral-moment for -s0, -weight(x)
  double GetSpectralMoment() {
    int N = getBinNumber();
    double sum = 0;
    for (int i=0; i <= N; i++) {
      sum += wRatio(s0_, weight_, i)*data_.GetSfm2(i);
    }
    return C::sTau/s0_/C::kBe*sum;
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
  std::function<double(double)> weight_;
  ExperimentalData data_;

  // get last included bin from s0_
  //
  // The experimental data is discrete. As we are summing over the data until
  // s0_ we need to know when to stop summing.
  int getBinNumber();

  std::vector<double> GetJacobianVector() {
    std::vector<double> jacobian(82);
    for (int i = 0; i < 80; i++) {
      int N = getBinNumber();
      if (i < N) {
        jacobian[i] = Constants::sTau/Constants::kBe/s0_
            *wRatio(s0_, weight_, i);
      } else {
        jacobian[i] = 0.;
      }
    }
    jacobian[80] = -GetSpectralMoment() / Constants::kBe;
    jacobian[81] = -GetSpectralMoment() / Constants::kRtauVex;
    return jacobian;
  }
};


#endif  // PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
