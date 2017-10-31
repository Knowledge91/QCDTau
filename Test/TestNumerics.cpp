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
  ASSERT_NEAR(Numerics::Integrate(realFunc, -0.2, 0.5), -0.602, maxError);
}

TEST(IntegralTest, Complex) {
  complex<double> a = -0.1+1i*0.1;
  complex<double> b = 0.3-1i*0.3;
  complex<double> r = Numerics::Integrate(complexFunc, a, b);

  ASSERT_NEAR(r.real(), -0.32, maxError);
  ASSERT_NEAR(r.imag(), -0.24, maxError);
}

TEST(MatrixInverterTest, Real) {
  //    1 2 3           -0.285714 0.261905 0.047619
  //A = 4 2 2 => A^-1 = 0.428571  0.190476 -0.238095
  //    5 1 7           0.142857  -0.214286 0.142857

  MatDoub_IO A(3, 3);
  A[0][0] = 1.; A[0][1] = 2.; A[0][2] = 3.;
  A[1][0] = 4.; A[1][1] = 2.; A[1][2] = 2.;
  A[2][0] = 5.; A[2][1] = 1.; A[2][2] = 7.;
  Numerics::gaussj(A);

  ASSERT_NEAR(A[0][0], -0.285714285714286, maxError);
  ASSERT_NEAR(A[0][1], 0.261904761904762, maxError);
  ASSERT_NEAR(A[0][2], 0.0476190476190476, maxError);
  ASSERT_NEAR(A[1][0], 0.428571428571429, maxError);
  ASSERT_NEAR(A[1][1], 0.190476190476190, maxError);
  ASSERT_NEAR(A[1][2], -0.238095238095238, maxError);
  ASSERT_NEAR(A[2][0], 0.142857142857143, maxError);
  ASSERT_NEAR(A[2][1], -0.214285714285714, maxError);
  ASSERT_NEAR(A[2][2], 0.142857142857143, maxError);
}
