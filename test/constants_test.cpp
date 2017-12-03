// Copyright 2017 Dirk Hornung

#include <cmath>
#include "gtest/gtest.h"
#include "../src/constants.h"

typedef Constants C;

class ConstantsTest : public ::testing::Test {
 protected:
  ConstantsTest() : constants_(C(3, 3, 4)) {}
  C constants_;
};


TEST_F(ConstantsTest, Zeta) {
  ASSERT_NEAR(constants_.getZeta(3), 1.20205690315959, C::maxError);
  ASSERT_NEAR(constants_.getZeta(5), 1.03692775514337, C::maxError);
  ASSERT_NEAR(constants_.getZeta(7), 1.00834927738192, C::maxError);
}

TEST_F(ConstantsTest, Beta) {
  ASSERT_NEAR(constants_.getBeta(1), 4.5, C::maxError);
  ASSERT_NEAR(constants_.getBeta(2), 8., C::maxError);
  ASSERT_NEAR(constants_.getBeta(3), 20.1197916666667, C::maxError*100);
  ASSERT_NEAR(constants_.getBeta(4), 94.4560791469040, C::maxError*1000);
}

TEST_F(ConstantsTest, Adler) {
  ASSERT_NEAR(constants_.getC(0, 1), 1., C::maxError);
  ASSERT_NEAR(constants_.getC(1, 1), 1., C::maxError);
  ASSERT_NEAR(constants_.getC(1, 2), 0., C::maxError);
  ASSERT_NEAR(constants_.getC(2, 1), 1.63982120489698, C::maxError*100);
  ASSERT_NEAR(constants_.getC(2, 2), -1.125, C::maxError);
  ASSERT_NEAR(constants_.getC(2, 3), 0., C::maxError);
  ASSERT_NEAR(constants_.getC(3, 1), 6.37101448310094, C::maxError*1000);
  ASSERT_NEAR(constants_.getC(3, 2), -5.6895977110182, C::maxError*1000);
  ASSERT_NEAR(constants_.getC(3, 3), 1.6875, C::maxError);
  ASSERT_NEAR(constants_.getC(3, 4), 0., C::maxError);
  ASSERT_NEAR(constants_.getC(4, 1), 49.0757000029480, C::maxError*100000);
  ASSERT_NEAR(constants_.getC(4, 2), -33.0914066167203, C::maxError*10000);
  ASSERT_NEAR(constants_.getC(4, 3), 15.8015948497910, C::maxError*1000);
  ASSERT_NEAR(constants_.getC(4, 4), -2.8476562500000, C::maxError);
  ASSERT_NEAR(constants_.getC(4, 5), 0., C::maxError);
  ASSERT_NEAR(constants_.getC(5, 1), 283., C::maxError);
  ASSERT_NEAR(constants_.getC(5, 2), -299.177187205152, C::maxError*100000);
  ASSERT_NEAR(constants_.getC(5, 3), 129.577532569234, C::maxError*10000);
  ASSERT_NEAR(constants_.getC(5, 4), -40.616088412030 , C::maxError*1000);
  ASSERT_NEAR(constants_.getC(5, 5), 5.1257812500000 , C::maxError);
}

TEST_F(ConstantsTest, various) {
  ASSERT_NEAR(C::sTau, 3.1570893124000001, C::maxError);
  ASSERT_NEAR(C::Be, 17.827, C::maxError);

  ASSERT_NEAR(C::Vud, 0.97425, C::maxError);
  ASSERT_NEAR(C::dVud, 0.00022, C::maxError);
  ASSERT_NEAR(C::fpi, 92.21e-3, C::maxError);
  ASSERT_NEAR(C::dfpi, 0.14e-3, C::maxError);
  ASSERT_NEAR(C::SEW, 1.0198, C::maxError);
  ASSERT_NEAR(C::dSEW, 0.0006, C::maxError);

  ASSERT_NEAR(C::pifac, 1.9494983309486447, C::maxError);
  ASSERT_NEAR(pow(C::dVud/C::Vud, 2.), 5.0992291959317588e-8, C::maxError);
  ASSERT_NEAR(pow(C::dSEW/C::SEW, 2.), 3.4615649558240843e-7, C::maxError);
  ASSERT_NEAR(4.*pow(C::dVud/C::Vud, 2.)+pow(C::dSEW/C::SEW, 2.)
              +4.*pow(C::dfpi/C::fpi, 2.), 9.7707434210522350e-6, C::maxError);
  ASSERT_NEAR(C::dpifac, 6.0937786116003652e-3, C::maxError);
}


