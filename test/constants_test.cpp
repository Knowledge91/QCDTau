// Copyright 2017 Dirk Hornung

#include <cmath>
#include "gtest/gtest.h"
#include "../src/constants.h"


class ConstantsTest : public ::testing::Test {
 protected:
  ConstantsTest() : constants_(Constants(3, 3, 4)) {}
  Constants constants_;
};


TEST_F(ConstantsTest, Zeta) {
  ASSERT_NEAR(constants_.getZeta(3), 1.20205690315959, Constants::maxError);
  ASSERT_NEAR(constants_.getZeta(5), 1.03692775514337, Constants::maxError);
  ASSERT_NEAR(constants_.getZeta(7), 1.00834927738192, Constants::maxError);
}

TEST_F(ConstantsTest, Beta) {
  ASSERT_NEAR(constants_.getBeta(1), 4.5, Constants::maxError);
  ASSERT_NEAR(constants_.getBeta(2), 8., Constants::maxError);
  ASSERT_NEAR(constants_.getBeta(3), 20.1197916666667, Constants::maxError*100);
  ASSERT_NEAR(constants_.getBeta(4), 94.4560791469040,
              Constants::maxError*1000);
}

TEST_F(ConstantsTest, Adler) {
  ASSERT_NEAR(constants_.getC(0, 1), 1., Constants::maxError);
  ASSERT_NEAR(constants_.getC(1, 1), 1., Constants::maxError);
  ASSERT_NEAR(constants_.getC(1, 2), 0., Constants::maxError);
  ASSERT_NEAR(constants_.getC(2, 1), 1.63982120489698, Constants::maxError*100);
  ASSERT_NEAR(constants_.getC(2, 2), -1.125, Constants::maxError);
  ASSERT_NEAR(constants_.getC(2, 3), 0., Constants::maxError);
  ASSERT_NEAR(constants_.getC(3, 1),
              6.37101448310094, Constants::maxError*1000);
  ASSERT_NEAR(constants_.getC(3, 2),
              -5.6895977110182, Constants::maxError*1000);
  ASSERT_NEAR(constants_.getC(3, 3), 1.6875, Constants::maxError);
  ASSERT_NEAR(constants_.getC(3, 4), 0., Constants::maxError);
  ASSERT_NEAR(constants_.getC(4, 1),
              49.0757000029480, Constants::maxError*100000);
  ASSERT_NEAR(constants_.getC(4, 2),
              -33.0914066167203, Constants::maxError*10000);
  ASSERT_NEAR(constants_.getC(4, 3),
              15.8015948497910, Constants::maxError*1000);
  ASSERT_NEAR(constants_.getC(4, 4), -2.8476562500000, Constants::maxError);
  ASSERT_NEAR(constants_.getC(4, 5), 0., Constants::maxError);
  ASSERT_NEAR(constants_.getC(5, 1), 283., Constants::maxError);
  ASSERT_NEAR(constants_.getC(5, 2),
              -299.177187205152, Constants::maxError*100000);
  ASSERT_NEAR(constants_.getC(5, 3),
              129.577532569234, Constants::maxError*10000);
  ASSERT_NEAR(constants_.getC(5, 4),
              -40.616088412030 , Constants::maxError*1000);
  ASSERT_NEAR(constants_.getC(5, 5), 5.1257812500000 , Constants::maxError);
}
