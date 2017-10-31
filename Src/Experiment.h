// Experiment.h
#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "Constants.h"
#include "Weights.h"
#include "Numerics.h"

#include <vector>

typedef Constants C;
typedef Weights W;
typedef Numerics N;

extern"C" {
  double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm, double *vpimhad, double *vprelepsp, double *vpreleptm, double *vpimlep, double *vpretopsp, double *vpretoptm, int *nrflag);
}
extern"C" {
  void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}

class Experiment {
 private:
  double sbin[80], dsbin[80], sfm2[80], derr[80], corerr[80][80];
  vector<vector<double>> errorMatrix;

  // get last included bin from s0
  int getBinNumber(double s0) {
    double pos = 0;
    for (int i=0; i<80; i++) {
      pos += dsbin[i];
      if(pos >= s0) {
        return i+1;
      }
    }
    return 80;
  }

  // \int_{ binStart }^{ binEnd} ds func(s)
  // Integrate weight functions over bin width
  template <typename Func>
  double binIntegral(int binNumber, Func &func) {
    double binStart = sbin[binNumber]-dsbin[binNumber]/2.;
    double binEnd = sbin[binNumber]-dsbin[binNumber]/2.;

    return N::Integrate(func, binStart, binEnd);
  }

  // fill Error Matrix
  void fillErrorMatrix() {
    //    cout << "derr " << derr[0] << " " << derr[1] << " " << derr[2] << endl;
    for (int i = 0; i < 80; i++) {
      for (int j = 0; j < 80; j++) {
        errorMatrix[i][j] = corerr[i][j]*derr[i]*derr[j]/100.;
      }
    }
  }

 public:
  Experiment()  : errorMatrix(80, vector<double>(80)) {
    // init Aleph Data
    aleph_vplusa_(sbin, dsbin, sfm2, derr, corerr);
    // init error Matrix
    fillErrorMatrix();
  }

  // get Spectral-moment for -s0, -weight(x)
  template <typename Func>
  double SpectralMoment(const double s0, Func &weight) {
    int N = getBinNumber(s0);
    cout << "s0s \t Nmax : " << s0 << "\t" << N << endl;
    cout << "s0/sTau \t" << s0/C::sTau << endl;
    cout << "sbin[0] \t" << sbin[0] << endl;
    cout << "B \t" << C::Be << endl;
    double sum = 0;
    for(int i=0; i<=N; i++) {

      //double wRatio = binIntegral(i, weight)/ binIntegral(i, W::WTau);
      double wRatio = s0/ C::sTau*(weight((sbin[i]-dsbin[i]/2.)/s0)-weight((sbin[i]+dsbin[i]/2.)/s0) )/ (W::WD00((sbin[i]-dsbin[i]/2.)/C::sTau) - W::WD00((sbin[i]+dsbin[i]/2.)/C::sTau));
      if (i==0) {
        cout << "weight \t" << weight((sbin[i]-dsbin[i]/2.)/s0)-weight((sbin[i]+dsbin[i]/2.)/s0)  << endl;
        cout << "weigthTau \t" << W::WD00((sbin[i]-dsbin[i]/2.)/C::sTau) - W::WD00((sbin[i]+dsbin[i]/2.)/C::sTau) << endl;
        cout << "wRatio \t" << wRatio << endl;
        cout << "factor \t" << C::sTau/s0/C::Be << endl;
        cout << "sfm2 \t" << sfm2[i] << endl;
      }

      sum += wRatio*sfm2[i];
    }
    return C::sTau/s0/C::Be*sum;
  }

  // output ErrorMatrix(i,j)
  void outErrorMatrix(int i, int j) {
    //cout << "corerr(2,1)" << corerr[0][0] << " " << corerr[0][1] << " " << corerr[1][1] << " " << corerr[1][0]  << endl;
    //cout << "ErrorMatrix(" << i << "," << j << ") = " << errorMatrix[i][j] << endl;
  }


};

#endif
