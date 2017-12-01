// Copyright 2017 Dirk Hornung

#include <iostream>
#include <cmath>
#include "CRunDec.h"
#include "Numerics.h"
#include "./constants.h"

using std::cout;
using std::endl;

typedef Constants C;

int main() {
  std::cout << std::setprecision(15);
  cout << C::pifac << endl;
  cout << M_PI << endl;
  cout << C::dpifac << endl;
  cout << 4.*std::pow(C::dVud/C::Vud, 2) + std::pow(C::dSEW/C::SEW, 2)+4.*std::pow(C::dfpi/C::fpi, 2) << endl;

  return(0);
}
