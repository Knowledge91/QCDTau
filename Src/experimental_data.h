// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_DATA_H_
#define PROGRAM_SRC_EXPERIMENTAL_DATA_H_

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>


// Access Fortran data provider
extern"C" {
  void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2,
                     double *derr, double (*corerr)[80]);
}
extern"C" {
  double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm,
                    double *vpimhad, double *vprelepsp, double *vpreleptm,
                    double *vpimlep, double *vpretopsp, double *vpretoptm,
                    int *nrflag);
}

namespace experimental_data_namespace {
  using std::vector;
  using boost::numeric::ublas::matrix;

class ExperimentalData {
 public:
  ExperimentalData() : sfm2_(80), sbin_(80), dsbin_(80), derr_(80),
      corerr_(80, 80), error_matrix_(80, 80) {
    double sbin[80], dsbin[80], sfm2[80], derr[80], corerr[80][80];
    aleph_vplusa_(sbin, dsbin, sfm2, derr, corerr);
    sfm2_.resize(80);
    for (int i = 0; i < GetNumberOfDataPoints(); i++) {
      sfm2_[i] = sfm2[i];
      sbin_[i] = sbin[i];
      dsbin_[i] = dsbin[i];
      // Normalize derr
      derr_[i] = 0.99363*derr[i];
      for (int j = 0; j < GetNumberOfDataPoints(); j++) {
        corerr_(i, j) = corerr[i][j];
        error_matrix_(i, j) = corerr[i][j]*derr_[i]*derr_[j]/100.;
       }
    }
  }

  int GetNumberOfDataPoints() const { return sfm2_.size(); }

  double GetSfm2(int i) const { return sfm2_[i]; }
  double GetSbin(int i) const { return sbin_[i]; }
  double GetDSbin(int i) const { return dsbin_[i]; }
  double GetDErr(int i) const { return derr_[i]; }
  double GetErrorMatrix(int i, int j) const { return error_matrix_(i, j); }
  matrix<double> GetCorrelationMatrix() const { return corerr_; }

 private:
  vector<double> sfm2_, sbin_, dsbin_, derr_;
  matrix<double> error_matrix_, corerr_;
};

}  // namespace  experimental_data_namespace

using experimental_data_namespace::ExperimentalData;

#endif  // PROGRAM_SRC_EXPERIMENTAL_DATA_H_
