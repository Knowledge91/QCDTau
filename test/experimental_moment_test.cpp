// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../src/experimental_moment.h"
#include "../src/weights.h"

typedef Weights W;

class ExperimentalMomentTest : public ::testing::Test {
 protected:
  ExperimentalMomentTest() :
      experimental_moment_(ExperimentalMoment(2.1000000000000001, W::WD00)) {}

  ExperimentalMoment experimental_moment_;
};

TEST_F(ExperimentalMomentTest, WRatio) {
  ASSERT_NEAR(experimental_moment_.wRatio(0),
              0.99930452590445074, C::maxError);
}

TEST_F(ExperimentalMomentTest, SpectralMoment) {
  ASSERT_NEAR(experimental_moment_.GetSpectralMoment(),
              2.2050182977666626e-5, C::maxError);
}

TEST_F(ExperimentalMomentTest, Jacobian) {
  ASSERT_NEAR(experimental_moment_.GetJacobianVector()[0],
                8.4272749429460669e-2, C::maxError);
  ASSERT_NEAR(experimental_moment_.GetJacobianVector()[80],
              -0.14706287667775386, C::maxError);
  ASSERT_NEAR(experimental_moment_.GetJacobianVector()[81],
              0.47026506524762296, C::maxError);
}

