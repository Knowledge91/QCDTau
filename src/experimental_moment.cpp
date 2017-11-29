// Copyright 2017 Dirk Hornung

#include "./experimental_moment.h"

using std::function;
using std::vector;

ExperimentalMoment::ExperimentalMoment(
    double s0, function<double(double)> weight) :
    s0_(s0), weight_(weight), data_(ExperimentalData()) {}

double ExperimentalMoment::GetSpectralMoment() {
  int N = getBinNumber();
  double sum = 0;
  for (int i=0; i <= N; i++) {
    sum += wRatio(i)*data_.GetSfm2(i);
  }
  return C::sTau/s0_/C::kBe*sum;
}

int ExperimentalMoment::getBinNumber() {
    double pos = 0;
    for (int i=0; i < 80; i++) {
      pos += data_.GetDSbin(i);
      if (pos >= s0_) {
        return i+1;
      }
    }
    return 80;
}

vector<double> ExperimentalMoment::GetJacobianVector() {
  vector<double> jacobian(82);
  for (int i = 0; i < 80; i++) {
    int N = getBinNumber();
    if (i < N) {
      jacobian[i] = Constants::sTau/Constants::kBe/s0_
          *wRatio(i);
    } else {
      jacobian[i] = 0.;
    }
  }
  jacobian[80] = -GetSpectralMoment() / Constants::kBe;
  jacobian[81] = -GetSpectralMoment() / Constants::kRtauVex;
  return jacobian;
}

double ExperimentalMoment::wRatio(int i) {
  return
      s0_/ C::sTau*
      (weight_((data_.GetSbin(i)-data_.GetDSbin(i)/2.)/s0_)
       -weight_((data_.GetSbin(i)+data_.GetDSbin(i)/2.)/s0_) )
      /
      (W::WD00((data_.GetSbin(i)-data_.GetDSbin(i)/2.)/C::sTau)
       - W::WD00((data_.GetSbin(i)+data_.GetDSbin(i)/2.)/C::sTau));
}

