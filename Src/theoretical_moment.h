// TheoreticalMoment.h
#ifndef THEORETICALMOMENT_H
#define THEORETICALMOMENT_H

#include "Constants.h"
#include "Numerics.h"
#include <cmath>
#include "CRunDec.h"
#include <complex>
#include "Numerics.h"


class TheoreticalMoment {
 public:
  TheoreticalMoment(Constants constants) : constants_(constants) {
    c_run_dec_ = new CRunDec(constants.nf);
  };

  complex<double> spectralMoment(double s0) {
    auto gamma = [s0] (complex<double> phi) {
      complex<double> i(0.0,1.0);
      return s0*exp(i*phi);
    };
    auto d0Gamma = [s0, this, gamma] (complex<double> x) {
      return D0(gamma(x));
    };

    complex<double> contourIntegral =  Numerics::Integrate(d0Gamma, complex<double>(0.,0.), complex<double>(2*M_PI,0.));
    return 3.*M_PI*contourIntegral;
  }

  complex<double> D0(complex<double> s) {
    complex<double> sum = 0;
    double amu = alphaMu(pow(Constants::kMu, 2));
    for (int n=1; n<=5; n++) {
      for (int k=1; k<=n; k++) {
        sum += k*pow(amu, n)*constants_.c[n][k]*pow(l(s, Constants::kMu), k-1);
      }
    }
    complex<double> d0 = constants_.nc/12./pow(M_PI, 2)*(constants_.c[0][1] + sum);
    if (first) {
      first = false;
      cout << endl;
      cout << "s \t" << s << endl;
      cout << "mu2 \t" << pow(Constants::kMu, 2) << endl;
      cout << "amu \t" << amu << endl;
      cout << "log(-s/mu2) \t" << l(s, Constants::kMu) << endl;
      cout << "D0 \t" << d0 << endl;
      cout << endl;
      for (int n=0; n<=5; n++) {
        for (int k=1; k<=n; k++) {
          cout << "c(" << n << "," << k << ") \t" << constants_.c[n][k] << endl;
        }
      }
    }
    return d0;
  }


 private:
  complex<double> l(complex<double> s, double mu) {
    complex<double> factor = -s/pow(mu,2);
    complex<double> result =  log(factor);
    //cout << "L\t" << result << endl;
    return result;
  }
  double alphaMu(double mu) {
    double alpha_s = c_run_dec_ -> AlphasExact(Constants::kAlphaSTau/M_PI, Constants::kSTau, mu, 3);
    cout << "alphaSTau, sTau" << Constants::kAlphaSTau/M_PI << " \t " << Constants::kSTau << endl;
    return alpha_s/M_PI;
  }

  bool first = true;
  Constants constants_;
  CRunDec* c_run_dec_;
};

#endif
