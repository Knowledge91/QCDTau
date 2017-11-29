// Copyright 2017 Dirk Hornung

#include <iostream>
#include <cmath>
#include "CRunDec.h"
#include "Numerics.h"
#include "./constants.h"


typedef Constants C;

int main() {
  int nf = 3;
  int nc = 3;
  int loops = 4;

  Constants constants(nc, nf, loops);

  std::cout << Constants::dpifac << std::endl;

  return(0);
}
