// Copyright 2017 Dirk Hornung

#include "./experimental_data.h"

ExperimentalData::ExperimentalData() : sfm2_(80), sbin_(80), dsbin_(80),
                                       derr_(80), corerr_(80, 80),
                                       error_matrix_(82, 82) {
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
  error_matrix_(80, 80) = pow(Constants::kDBe, 2);
  error_matrix_(81, 81) = 1.;
}

void ExperimentalData::test() {
  std::cout << Constants::dpifac << std::endl;
}
