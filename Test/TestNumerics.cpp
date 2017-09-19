#include <complex>
#include "gtest/gtest.h"
#include "Numerics.h"

using namespace std;
using namespace std::complex_literals;

double realFunc(double x) {
  return 3*x*x-7*x;
}
complex<double> complexFunc(complex<double> x) {
  return 3.*x-1i*4.*x;
}


const double maxError = 1e-15;

TEST(IntegralTest, Real) {
  ASSERT_NEAR(Numerics::qgauss(realFunc, -0.2, 0.5), -0.602, maxError);
  complex<double> a = -0.1+1i*0.1;
  complex<double> b = 0.3-1i*0.3;
  complex<double> r = Numerics::qgauss(complexFunc, a, b);

  ASSERT_NEAR(r.real(), -0.32, maxError);
  ASSERT_NEAR(r.imag(), -0.24, maxError);
}
