// Copyright 2017 Dirk Hornung

#include <vector>
#include "gtest/gtest.h"
#include "../src/constants.h"
#include "../src/chisquared.h"

using std::vector;
class ChisquaredTest : public ::testing::Test {
 protected:
  ChisquaredTest() :
      chisquared_(Chisquared(Constants(3, 3, 4),
                 {C::sTau, 3., 2.8, 2.6, 2.4, 2.3, 2.2, 2.1, 2.})) {}

  Chisquared chisquared_;
};


TEST_F(ChisquaredTest, JacobianMatrix) {
  ASSERT_NEAR(chisquared_.GetJacobianMatrix(0, 0), 5.60946878330622e-2,
              C::maxError*1e1);
  ASSERT_DOUBLE_EQ(chisquared_.GetJacobianMatrix(0, 7),
              8.4272749429460669e-2);
  ASSERT_NEAR(chisquared_.GetJacobianMatrix(21, 4), 6.856547153793633e-2,
              C::maxError*1e1);
}

TEST_F(ChisquaredTest, CovarianceMatrix) {
  ASSERT_NEAR(chisquared_.GetCovarianceMatrix(0, 0),
              9.2174853727193312e-5, C::maxError*1e10);
  // ASSERT_NEAR(chisquared_.GetCovarianceMatrix(1, 1),
  //              6.6522888710431944e-5, C::maxError);
}
