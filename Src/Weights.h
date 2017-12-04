// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_WEIGHTS_H_
#define PROGRAM_SRC_WEIGHTS_H_

#include "./constants.h"

typedef Constants C;

class Weights {
 public:
  static double WTau(double x) {
    return (1.-x)*(1.-x)*(1.+2.*x);
  }
  static double WD00(double x) {
    return pow(1.-x, 3)*(1.+x);
  }
  static double WR00(double x) {
    return pow(1.-x, 2)*(1.+2.*x);
  }
};

#endif  // PROGRAM_SRC_WEIGHTS_H_
