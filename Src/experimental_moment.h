// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
#define PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_

#include <vector>
#include "Constants.h"
#include "Weights.h"
#include "./Numerics.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>



typedef Constants C;
typedef Weights W;
typedef Numerics N;

extern"C" {
  double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm, double *vpimhad, double *vprelepsp, double *vpreleptm, double *vpimlep, double *vpretopsp, double *vpretoptm, int *nrflag);
}
extern"C" {
  void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}

namespace Ublas {

  using boost::numeric::ublas::matrix;
  using std::vector;
  using std::cout;

class ExperimentalMoment {
 public:
  ExperimentalMoment() : errorMatrix(80, 80), covarianceMatrix(80, 80),
                         jacobian(80) {
      // init Aleph Data
      aleph_vplusa_(sbin, dsbin, sfm2, derr, corerr);
      // init covariant matrix
      fillErrorMatrix();
      // Numerics::outputMatrix(errorMatrix);
      fillJacobian(1., W::WTau);
      Numerics::outputVector(jacobian);
      fillCovarianceMatrix();
      Numerics::outputMatrix(covarianceMatrix);
  }

  // get Spectral-moment for -s0, -weight(x)
  double SpectralMoment(const double s0, function<double(double)> weight) {
    int N = getBinNumber(s0);
    cout << "s0s \t Nmax : " << s0 << "\t" << N << endl;
    cout << "s0/sTau \t" << s0/C::sTau << endl;
    cout << "sbin[0] \t" << sbin[0] << endl;
    cout << "B \t" << C::Be << endl;
    double sum = 0;
    for(int i=0; i<=N; i++) {

      //double wRatio = binIntegral(i, weight)/ binIntegral(i, W::WTau);
      // double wRatio = s0/ C::sTau*(weight((sbin[i]-dsbin[i]/2.)/s0)-weight((sbin[i]+dsbin[i]/2.)/s0) )/ (W::WD00((sbin[i]-dsbin[i]/2.)/C::sTau) - W::WD00((sbin[i]+dsbin[i]/2.)/C::sTau));
      if (i==0) {
        cout << "weight \t" << weight((sbin[i]-dsbin[i]/2.)/s0)-weight((sbin[i]+dsbin[i]/2.)/s0)  << endl;
        cout << "weigthTau \t" << W::WD00((sbin[i]-dsbin[i]/2.)/C::sTau) - W::WD00((sbin[i]+dsbin[i]/2.)/C::sTau) << endl;
	//        cout << "wRatio \t" << wRatio << endl;
        cout << "factor \t" << C::sTau/s0/C::Be << endl;
        cout << "sfm2 \t" << sfm2[i] << endl;
      }

      sum += wRatio(s0, weight, i)*sfm2[i];
    }
    return C::sTau/s0/C::Be*sum;
  }

 private:
  double sbin[80], dsbin[80], sfm2[80], derr[80], corerr[80][80];
  matrix<double> errorMatrix, covarianceMatrix;
  vector<double> jacobian;

  // get last included bin from s0
  int getBinNumber(double s0) {
    double pos = 0;
    for (int i=0; i < 80; i++) {
      pos += dsbin[i];
      if (pos >= s0) {
        return i+1;
      }
    }
    return 80;
  }

  // fill Error Matrix
  void fillErrorMatrix() {
    for (int i = 0; i < 80; i++) {
      for (int j = 0; j < 80; j++) {
        errorMatrix(i, j) = corerr[i][j]*derr[i]*derr[j]/100.;
      }
    }
  }

  void fillJacobian(double s0, function<double(double)> weight) {
    for (int i = 0; i < 80; i++) {
      int N = getBinNumber(s0);
      if (i < N) {
        jacobian[i] = Constants::sTau/Constants::Be/s0 * wRatio(s0, weight, i);
      } else {
        jacobian[i] = 0.;
      }
    }
  }

  void fillCovarianceMatrix() {
    for (int i = 0; i < 80; i++) {
      for (int j = 0; j < 80; j++) {
        covarianceMatrix(i, j) = jacobian[i] * errorMatrix(i, j) * jacobian[j];
      }
    }
  }

  template <typename Func>
  double wRatio(double s0, Func weight, int i) {
    return s0/ C::sTau*(weight((sbin[i]-dsbin[i]/2.)/s0)-weight((sbin[i]+dsbin[i]/2.)/s0) )/ (W::WD00((sbin[i]-dsbin[i]/2.)/C::sTau) - W::WD00((sbin[i]+dsbin[i]/2.)/C::sTau));
  }
};

}  // namespace Ublas

using Ublas::ExperimentalMoment;

#endif  // PROGRAM_SRC_EXPERIMENTAL_MOMENT_H_
