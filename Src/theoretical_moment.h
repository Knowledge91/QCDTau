// TheoreticalMoment.h
#ifndef THEORETICALMOMENT_H
#define THEORETICALMOMENT_H

#include "./constants.h"
#include "./numerics.h"
#include <cmath>
#include "CRunDec.h"
#include <complex>


class TheoreticalMoment {
 public:
  TheoreticalMoment(Constants constants) : constants_(constants) {
    c_run_dec_ = new CRunDec(constants.getNf());
  };

  std::complex<double> spectralMoment(double s0) {
    auto gamma = [s0] (std::complex<double> phi) {
      std::complex<double> i(0.0,1.0);
      return s0*exp(i*phi);
    };
    auto d0Gamma = [s0, this, gamma] (std::complex<double> x) {
      return D0(gamma(x));
    };

    std::complex<double> contourIntegral =  Numerics::Integrate(d0Gamma, std::complex<double>(0.,0.), std::complex<double>(2*M_PI,0.));
    return 3.*M_PI*contourIntegral;
  }

  std::complex<double> D0(std::complex<double> s) {
    std::complex<double> sum = 0;
    double amu = alphaMu(pow(Constants::kMu, 2));
    for (int n=1; n<=5; n++) {
      for (int k=1; k<=n; k++) {
        sum += k*pow(amu, n)*constants_.getC(n, k)*pow(l(s, Constants::kMu), k-1);
      }
    }
    std::complex<double> d0 = constants_.getNc()/12./pow(M_PI, 2)*(constants_.getC(0, 1) + sum);
    if (first) {
      first = false;
      std::cout << std::endl;
      std::cout << "s \t" << s << std::endl;
      std::cout << "mu2 \t" << pow(Constants::kMu, 2) << std::endl;
      std::cout << "amu \t" << amu << std::endl;
      std::cout << "log(-s/mu2) \t" << l(s, Constants::kMu) << std::endl;
      std::cout << "D0 \t" << d0 << std::endl;
      std::cout << std::endl;
      for (int n=0; n<=5; n++) {
        for (int k=1; k<=n; k++) {
          std::cout << "c(" << n << "," << k << ") \t" << constants_.getC(n, k) <<
              std::endl;
        }
      }
    }
    return d0;
  }


 private:
  std::complex<double> l(std::complex<double> s, double mu) {
    std::complex<double> factor = -s/pow(mu,2);
    std::complex<double> result =  log(factor);
    //cout << "L\t" << result << endl;
    return result;
  }
  double alphaMu(double mu) {
    double alpha_s = c_run_dec_ -> AlphasExact(Constants::kAlphaSTau/M_PI, Constants::kSTau, mu, 3);
    std::cout << "alphaSTau, sTau" << Constants::kAlphaSTau/M_PI << " \t " << Constants::kSTau << std::endl;
    return alpha_s/M_PI;
  }

  bool first = true;
  Constants constants_;
  CRunDec* c_run_dec_;
};

#endif
