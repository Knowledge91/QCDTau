// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_DATA_H_
#define PROGRAM_SRC_EXPERIMENTAL_DATA_H_

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "./constants.h"


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

class ExperimentalData {
 public:
  ExperimentalData();

  int GetNumberOfDataPoints() const { return sfm2_.size(); }

  double GetSfm2(int i) const { return sfm2_[i]; }
  double GetSbin(int i) const { return sbin_[i]; }
  double GetDSbin(int i) const { return dsbin_[i]; }
  double GetDErr(int i) const { return derr_[i]; }
  double GetErrorMatrix(int i, int j) const { return error_matrix_(i, j); }
  boost::numeric::ublas::matrix<double> GetCorrelationMatrix() const {
    return corerr_;
  }

 private:
  std::vector<double> sfm2_, sbin_, dsbin_, derr_;
  boost::numeric::ublas::matrix<double> error_matrix_, corerr_;
};



#endif  // PROGRAM_SRC_EXPERIMENTAL_DATA_H_
