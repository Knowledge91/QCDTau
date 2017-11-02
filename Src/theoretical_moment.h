// TheoreticalMoment.h
#ifndef THEORETICALMOMENT_H
#define THEORETICALMOMENT_H

#include "Constants.h"
#include "Numerics.h"
#include <cmath>
#include "CRunDec.h"
#include <complex>


class TheoreticalMoment {
 public:
  TheoreticalMoment(Constants constants) : constants_(constants) {
    c_run_dec_ = new CRunDec(constants.nf);
    cout << "Init TheoreticalMoments \t" << alphaMu(Constants::sTau) << endl;
  };

  complex<double> D0(double s) {
    complex<double> sum = 0;
    for (int n=1; n<=5; n++) {
      for (int k=1; k<=n; k++) {
        sum += k*pow(alphaMu(s), n)*constants_.c[n][k]*pow(l(s, constants_.kMu), k-1);
      }
    }
    return constants_.nc/12./pow(M_PI, 2)*(constants_.c[0][1] + sum);
  }
 private:
  complex<double> l(double s, double mu) {
    complex<double> factor = -s/pow(mu,2);
    complex<double> result =  log(factor);
    //cout << "L\t" << result << endl;
    return result;
  }
  double alphaMu(double s) {
    double alpha_s = c_run_dec_ -> AlphasExact(constants_.kAlphaSMz, constants_.kMz, s, 4);
    //cout << "amu(" << s << ")" << alpha_s/M_PI << endl;
    return alpha_s/M_PI;
  }

  Constants constants_;
  CRunDec* c_run_dec_;
};

#endif
