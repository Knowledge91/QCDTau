// Copyright 2017 Dirk Hornung

#include "./experimental_moment.h"

using std::function;
using std::vector;
using std::cout;
using std::endl;

ExperimentalMoment::ExperimentalMoment(
    double s0, function<double(double)> weight) :
    s0_(s0), weight_(weight), data_(ExperimentalData()) {}

double ExperimentalMoment::GetSpectralMoment() {
  double sum = 0;
  for (int i=0; i <= GetBinNumber(); i++) {
    sum += C::sTau/s0_/C::Be*wRatio(i)*data_.GetScaledSfm2(i);
  }
  return sum + GetPiMoment();
}

double ExperimentalMoment::GetPiMoment() {
  // Axialvector pion-pole contribution
  double momA = C::pifac/s0_*W::WR00(pow(C::mpim, 2)/s0_);
  // Pseudoscalar pion-pole contribution
  double momP = momA*(-2.*pow(C::mpim, 2)/(C::sTau+2.*pow(C::mpim, 2)));
  return momA + momP;
}

int ExperimentalMoment::GetBinNumber() {
    double pos = 0;
    for (int i=0; i < data_.GetNumberOfDataPoints(); i++) {
      if (fabs(s0_-6.e-6-data_.GetSbin(i)) < data_.GetDSbin(i)/2.) {
        return i;
      }
    }
    return 80;
}

vector<double> ExperimentalMoment::GetJacobianVector() {
  vector<double> jacobian(82);
  for (int i = 0; i < 80; i++) {
    int N = GetBinNumber();
    if (i < N) {
      jacobian[i] = Constants::sTau/Constants::Be/s0_
          *wRatio(i);
    } else {
      jacobian[i] = 0.;
    }
  }
  jacobian[80] = (GetPiMoment()-GetSpectralMoment()) / Constants::Be;
  jacobian[81] = GetPiMoment() / Constants::pifac;
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

