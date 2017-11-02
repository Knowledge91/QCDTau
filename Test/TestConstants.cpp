#include "gtest/gtest.h"
#include <cmath>
#include "Constants.h"


// Tests for nf = 3
Constants C(3,3,4);

TEST(ConstantsTest, Zeta) {
  ASSERT_NEAR(C.zeta[3], 1.20205690315959, Constants::maxError);
  ASSERT_NEAR(C.zeta[5], 1.03692775514337, Constants::maxError);
  ASSERT_NEAR(C.zeta[7], 1.00834927738192, Constants::maxError);
}

TEST(ConstantsTest, Beta) {
  ASSERT_NEAR(C.beta[1], 4.5, Constants::maxError);
  ASSERT_NEAR(C.beta[2], 8., Constants::maxError);
  ASSERT_NEAR(C.beta[3], 20.1197916666667, Constants::maxError*100);
  ASSERT_NEAR(C.beta[4], 94.4560791469040, Constants::maxError*1000);
}

TEST(ConstantsTest, Adler) {
  ASSERT_NEAR(C.c[0][1], 1., Constants::maxError);
  ASSERT_NEAR(C.c[1][1], 1., Constants::maxError);
  ASSERT_NEAR(C.c[1][2], 0., Constants::maxError);
  ASSERT_NEAR(C.c[2][1], 1.63982120489698, Constants::maxError*100);
  ASSERT_NEAR(C.c[2][2], 1.125, Constants::maxError);
  ASSERT_NEAR(C.c[2][3], 0., Constants::maxError);
  ASSERT_NEAR(C.c[3][1], 6.37101448310094, Constants::maxError*1000);
  ASSERT_NEAR(C.c[3][2], -5.6895977110182, Constants::maxError*1000);
  ASSERT_NEAR(C.c[3][3], 1.6875, Constants::maxError);
  ASSERT_NEAR(C.c[3][4], 0., Constants::maxError);
  ASSERT_NEAR(C.c[4][1], 49.0757000029480, Constants::maxError*1000);
  ASSERT_NEAR(C.c[4][2], -33.0914066167203, Constants::maxError*10000);
  ASSERT_NEAR(C.c[4][3], 15.8015948497910, Constants::maxError*1000);
  ASSERT_NEAR(C.c[4][4], -2.8476562500000, Constants::maxError);
  ASSERT_NEAR(C.c[4][5], 0., Constants::maxError);


}
