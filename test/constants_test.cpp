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
  // compared to Matthias: num_const.f90
  ASSERT_DOUBLE_EQ(constants_.GetZeta(3), 1.202056903159594285);
  ASSERT_DOUBLE_EQ(constants_.GetZeta(5), 1.036927755143369926);
  ASSERT_DOUBLE_EQ(constants_.GetZeta(7), 1.008349277381922827);
}

TEST_F(ConstantsTest, Beta) {
  ASSERT_DOUBLE_EQ(constants_.GetBeta(1), 4.5);
  ASSERT_DOUBLE_EQ(constants_.GetBeta(2), 8.);
  ASSERT_DOUBLE_EQ(constants_.GetBeta(3), 20.119791666666668);
  ASSERT_DOUBLE_EQ(constants_.GetBeta(4), 94.456079146903988);
}


TEST_F(ConstantsTest, Adler) {
  // compare to coefficients.nb
  ASSERT_DOUBLE_EQ(constants_.GetC(0, 1), 1.);
  ASSERT_DOUBLE_EQ(constants_.GetC(1, 1), 1.);
  ASSERT_DOUBLE_EQ(constants_.GetC(1, 2), 0.);
  ASSERT_DOUBLE_EQ(constants_.GetC(2, 1), 1.639821204896986);
  ASSERT_DOUBLE_EQ(constants_.GetC(2, 2), -1.125);
  ASSERT_DOUBLE_EQ(constants_.GetC(2, 3), 0.);
  ASSERT_DOUBLE_EQ(constants_.GetC(3, 1), 6.371014483100957);
  ASSERT_DOUBLE_EQ(constants_.GetC(3, 2), -5.68959771101822);
  ASSERT_DOUBLE_EQ(constants_.GetC(3, 3), 1.6875);
  ASSERT_DOUBLE_EQ(constants_.GetC(3, 4), 0.);
  ASSERT_DOUBLE_EQ(constants_.GetC(4, 1), 49.075700002948054);
  ASSERT_DOUBLE_EQ(constants_.GetC(4, 2), -33.091406616720334);
  ASSERT_DOUBLE_EQ(constants_.GetC(4, 3), 15.8015948497910);
  ASSERT_DOUBLE_EQ(constants_.GetC(4, 4), -2.8476562500000);
  ASSERT_DOUBLE_EQ(constants_.GetC(4, 5), 0.);
  ASSERT_DOUBLE_EQ(constants_.GetC(5, 1), 283.);
  ASSERT_DOUBLE_EQ(constants_.GetC(5, 2), -299.17718720515285);
  ASSERT_DOUBLE_EQ(constants_.GetC(5, 3), 129.577532569234);
  ASSERT_DOUBLE_EQ(constants_.GetC(5, 4), -40.61608841202971);
  ASSERT_DOUBLE_EQ(constants_.GetC(5, 5), 5.1257812500000);

  // compare to Matthias:
  EXPECT_DOUBLE_EQ(constants_.GetC(0, 1), 1.);
  EXPECT_DOUBLE_EQ(constants_.GetC(1, 1), 1.);
  // EXPECT_DOUBLE_EQ(constants_.GetC(2, 1), 1.63982120489698476474);
  //  EXPECT_DOUBLE_EQ(constants_.GetC(3, 1), 6.37101448310094071138);
  //  EXPECT_DOUBLE_EQ(constants_.GetC(4, 1), 49.07570000294798513221);
  EXPECT_DOUBLE_EQ(constants_.GetC(5, 1), 283.);
}

TEST_F(ConstantsTest, various) {
  ASSERT_NEAR(C::sTau, 3.1570893124000001, C::maxError);
  ASSERT_NEAR(C::Be, 17.827000000000002, C::maxError);

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


