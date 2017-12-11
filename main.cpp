// Copyright 2017 Dirk Hornung

#include <iostream>
#include <typeinfo>
#include <cmath>
#include "CRunDec.h"
#include "src/numerics.h"
#include "src/constants.h"
#include "src/weights.h"
#include "src/experimental_moment.h"



using std::cout;
using std::endl;
typedef Constants C;
typedef Weights W;


int main() {
  cout << std::setprecision(17);
  C constants(3, 3, 3);

  //  cout << typeid(pow(3, 5)).name() << endl;
  ExperimentalMoment experimental_moment(2.1, W::WD00);
  //  cout << C::sTau/2.1/C::Be - 8.4331399733416662e-2 << endl;

  cout << 149753./768. - 1103./4.*constants.GetZeta(3) + 275./6.*constants.GetZeta(5) + 88.95045018552034  << endl;
  cout << constants.GetZeta(5) - 1.03692775514337 << endl;
  cout << 275./6.*constants.GetZeta(5)  - 47.52585544407113 << endl;

  return(0);
}
