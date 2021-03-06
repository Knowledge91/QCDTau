// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../src/experimental_moment.h"
#include "../src/weights.h"

typedef Weights W;

class ExperimentalMomentTest : public ::testing::Test {
 protected:
  ExperimentalMomentTest() :
      experimental_moment_(ExperimentalMoment(2.1, W::WD00)),
      experimental_moment2_(ExperimentalMoment(2., W::WD00)) {}

  ExperimentalMoment experimental_moment_;     // matthias s0s(1,8)
  ExperimentalMoment experimental_moment2_;    // matthias s0s(1,9)
};

TEST_F(ExperimentalMomentTest, WRatio) {
  ASSERT_DOUBLE_EQ(experimental_moment_.wRatio(0), 0.99930452590445074);
  ASSERT_NEAR(experimental_moment2_.wRatio(5), 0.98125247104896496,
              C::maxError*1e1);
}

TEST_F(ExperimentalMomentTest, BinNumber) {
  ASSERT_EQ(experimental_moment_.GetBinNumber(), 71);
  ASSERT_EQ(experimental_moment2_.GetBinNumber(), 70);
}

TEST_F(ExperimentalMomentTest, PiMoment) {
  ASSERT_NEAR(experimental_moment_.GetPiMoment(),
              0.91678095980369645, C::maxError);
}

TEST_F(ExperimentalMomentTest, SpectralMoment) {
  // floating point summation introduced e-1 error
  ASSERT_NEAR(experimental_moment_.GetSpectralMoment(),
              3.5421269656941181, C::maxError*1e1);
}

TEST_F(ExperimentalMomentTest, Jacobian) {
  ASSERT_DOUBLE_EQ(Constants::sTau/(Constants::Be*2.1),
                   Constants::sTau/Constants::Be/2.1);
  ASSERT_DOUBLE_EQ(experimental_moment_.GetJacobianVector()[0],
                8.4272749429460669e-2);
  // jacobian[80] depends on spectralMoment
  ASSERT_NEAR(experimental_moment_.GetJacobianVector()[80],
              -0.14726796465419989, C::maxError*1e1);
  ASSERT_DOUBLE_EQ(experimental_moment_.GetJacobianVector()[81],
              0.47026506524762296);
}

