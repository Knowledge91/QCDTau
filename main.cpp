// Copyright 2017 Dirk Hornung

#include <iostream>
#include <cmath>
#include "CRunDec.h"
#include "src/numerics.h"
#include "src/constants.h"
#include "src/weights.h"
#include "src/experimental_moment.h"

typedef Constants C;
typedef Weights W;



using std::cout;
using std::endl;

int main() {
  std::cout << std::setprecision(15);

  ExperimentalMoment experimental_moment(2.1, W::WD00);
  cout << experimental_moment.GetSpectralMoment() - 47.033957394180902 << endl;
  cout << C::sTau/2.1/C::Be - 8.4331399733416662e-2 << endl;

  return(0);
}
