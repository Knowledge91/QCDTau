// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_EXPERIMENTAL_DATA_H_
#define PROGRAM_SRC_EXPERIMENTAL_DATA_H_

#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

extern"C" {
  void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}

namespace experimental_data_namespace {
  using std::vector;
  using boost::numeric::ublas::matrix;

class ExperimentalData {
 public:
ExperimentalData() : sfm2_(80), sbin_(80), dsbin_(80), derr_(80),
      corerr_(80, 80) {
    n_data_points_ = 80;
    double sbin[80], dsbin[80], sfm2[80], derr[80], corerr[80][80];
    aleph_vplusa_(sbin, dsbin, sfm2, derr, corerr);
    sfm2_.resize(80);
    for (int i = 0; i < n_data_points_; i++) {
      sfm2_[i] = sfm2[i];
      sbin_[i] = sbin[i];
      dsbin_[i] = dsbin[i];
      derr_[i] = derr[i];
      for (int j = 0; j < n_data_points_; j++) {
        corerr_(i, j) = corerr[i][j];
      }
    }

    std::cout << sfm2[1] << std::endl;
  }

 private:
  int n_data_points_;
  vector<double> sfm2_, sbin_, dsbin_, derr_;
  matrix<double> corerr_;
};

}  // namespace  experimental_data_namespace

using experimental_data_namespace::ExperimentalData;

#endif  // PROGRAM_SRC_EXPERIMENTAL_DATA_H_
