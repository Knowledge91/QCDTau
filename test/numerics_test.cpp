// Copyright 2017 Dirk Hornung

#include <complex>
#include "gtest/gtest.h"
#include "../src/numerics.h"

using std::complex;
using std::complex_literals::operator""i;
using boost::numeric::ublas::matrix;

typedef Constants C;

double realFunc(double x) {
  return 3*x*x-7*x;
}
complex<double> complexFunc(complex<double> x) {
  return 3.*x-1i*4.*x;
}


TEST(NumericsTest, RealIntegral) {
  ASSERT_NEAR(Numerics::Integrate(realFunc, -0.2, 0.5), -0.602, C::maxError);
}

TEST(NumericsTest, ComplexIntegral) {
  complex<double> a = -0.1+1i*0.1;
  complex<double> b = 0.3-1i*0.3;
  complex<double> r = Numerics::Integrate(complexFunc, a, b);

  ASSERT_NEAR(r.real(), -0.32, C::maxError);
  ASSERT_NEAR(r.imag(), -0.24, C::maxError);
}

TEST(NumericsTest, MatrixInverse) {
  //      1 2 3           -0.285714 0.261905 0.047619
  //  A = 4 2 2 => A^-1 = 0.428571  0.190476 -0.238095
  //      5 1 7           0.142857  -0.214286 0.142857

  matrix<double>  A(3, 3);
  matrix<double> AInv(3, 3);
  A(0, 0) = 1.; A(0, 1) = 2.; A(0, 2) = 3.;
  A(1, 0) = 4.; A(1, 1) = 2.; A(1, 2) = 2.;
  A(2, 0) = 5.; A(2, 1) = 1.; A(2, 2) = 7.;
  Numerics::InvertMatrix(A, AInv);

  ASSERT_NEAR(AInv(0, 0), -0.285714285714286, C::maxError);
  ASSERT_NEAR(AInv(0, 1), 0.261904761904762, C::maxError);
  ASSERT_NEAR(AInv(0, 2), 0.0476190476190476, C::maxError);
  ASSERT_NEAR(AInv(1, 0), 0.428571428571429, C::maxError);
  ASSERT_NEAR(AInv(1, 1), 0.190476190476190, C::maxError);
  ASSERT_NEAR(AInv(1, 2), -0.238095238095238, C::maxError);
  ASSERT_NEAR(AInv(2, 0), 0.142857142857143, C::maxError);
  ASSERT_NEAR(AInv(2, 1), -0.214285714285714, C::maxError);
  ASSERT_NEAR(AInv(2, 2), 0.142857142857143, C::maxError);

  matrix<double> B(3, 3);
  matrix<double> BInv(3, 3);
  B(0, 0) = 1.e-6; B(0, 1) = 4.e-6; B(0, 2) = 5.e-3;
  B(1, 0) = 2.e-7; B(1, 1) = 2.e-8; B(1, 2) = 1.e-6;
  B(2, 0) = 3.e-6; B(2, 1) = 2.e-6; B(2, 2) = 7.e-6;
  Numerics::InvertMatrix(B, BInv);

  double inv_err = C::maxError*1e7;
  ASSERT_NEAR(BInv(0, 0), -1091.20349185117, inv_err);
  ASSERT_NEAR(BInv(0, 1), 5.85025872082791e6, inv_err);
  ASSERT_NEAR(BInv(0, 2), -56320.1802245767, inv_err);
  ASSERT_NEAR(BInv(1, 0), 938.669670409612, inv_err);
  ASSERT_NEAR(BInv(1, 1), -8.79592148028207e6, inv_err);
  ASSERT_NEAR(BInv(1, 2), 586081.875462001, inv_err);
  ASSERT_NEAR(BInv(2, 0), 199.467304962043, inv_err);
  ASSERT_NEAR(BInv(2, 1), 5866.68544006007, inv_err);
  ASSERT_NEAR(BInv(2, 2), -457.601464324686, inv_err);
}
