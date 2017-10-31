// Weights.h
#ifndef WEIGHTS_H
#define WEIGHTS_H

#include "Constants.h"

typedef Constants C;

class Weights {
 public:
  static double WTau(double x) {
    return (1-x)*(1-x)*(1+2*x);
  }
  static double WD00(double x) {
    return pow(1-x, 3)*(1+x);
  }
};

#endif
