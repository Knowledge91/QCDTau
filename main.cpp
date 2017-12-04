// Copyright 2017 Dirk Hornung

#include <iostream>
#include <cmath>
#include "CRunDec.h"
#include "Numerics.h"
#include "src/constants.h"
#include "src/weights.h"
#include "src/experimental_moment.h"

typedef Weights W;

using std::cout;
using std::endl;

typedef Constants C;

int main() {
  std::cout << std::setprecision(15);

  ExperimentalMoment experimental_moment(2.1, W::WD00);
  cout << experimental_moment.GetSpectralMoment() - 47.033957394180902 << endl;
  cout << C::sTau/2.1/C::Be - 8.4331399733416662e-2 << endl;

  return(0);
}
