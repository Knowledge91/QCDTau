// Copyright 2017 Dirk Hornung

#include "./theoretical_moment.h"

using std::complex;
typedef complex<double> dComplex;

dComplex TheoreticalMoment::GetSpectralMoment(double s0) {
  auto gamma = [s0] (dComplex phi) {
    dComplex i(0.0, 1.0);
      return s0*exp(i*phi);
    };
  auto d0Gamma = [s0, this, gamma] (dComplex x) {
    return GetD0(gamma(x));
  };

  dComplex contourIntegral =
      Numerics::Integrate(d0Gamma, dComplex(0., 0.), dComplex(2*M_PI, 0.));

  return 3.*M_PI*contourIntegral;
}

dComplex TheoreticalMoment::GetD0(dComplex s) {
  dComplex sum = 0;
  double amu = GetAlphaMu(pow(Constants::kMu, 2));
  for (int n=1; n <= 5; n++) {
    for (int k=1; k <= n; k++) {
      sum += k*pow(amu, n)*constants_.getC(n, k)
          *pow(GetL(s, Constants::kMu), k-1);
    }
  }
  dComplex d0 = constants_.getNc()/12./pow(M_PI, 2)
      *(constants_.getC(0, 1) + sum);
  for (int n=0; n <= 5; n++) {
    for (int k=1; k <= n; k++) {
      std::cout << "c(" << n << "," << k << ") \t"
                << constants_.getC(n, k) << std::endl;
    }
  }
  return d0;
}

dComplex TheoreticalMoment::GetL(dComplex s, double mu) {
  dComplex factor = -s/pow(mu, 2);
  dComplex result =  log(factor);
  return result;
}

double TheoreticalMoment::GetAlphaMu(double mu)  {
  double alpha_s = c_run_dec_ -> AlphasExact(Constants::kAlphaSTau/M_PI,
                                             Constants::kSTau, mu, 3);
  return alpha_s/M_PI;
}


