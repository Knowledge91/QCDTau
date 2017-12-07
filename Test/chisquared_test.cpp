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
  ASSERT_NEAR(chisquared_.GetJacobianMatrix(62, 2), 5.01150477053736561e-2,
              C::maxError*1e1);
}

TEST_F(ChisquaredTest, CovarianceMatrix) {
  double test1 = chisquared_.jacobian_matrix_(22, 7)*
      chisquared_.data_.GetErrorMatrix(1, 3)*chisquared_.jacobian_matrix_(1, 2);
  double test2 = chisquared_.jacobian_matrix_(48, 3)*
      chisquared_.data_.GetErrorMatrix(3, 5)*chisquared_.jacobian_matrix_(8, 8);
  ASSERT_NEAR(test1, -7.9851491495089808e-10, C::maxError);
  ASSERT_NEAR(test2, 1.3672464839786044e-8, C::maxError);
  // !!! Covariance matrix is only to 1e4 accurate!
  ASSERT_NEAR(chisquared_.GetCovarianceMatrix(0, 0),
              1.3667648017148091e-4, C::maxError*1e10);
  ASSERT_NEAR(chisquared_.GetCovarianceMatrix(1, 1),
              6.6522888710431944e-5, C::maxError*1e11);
}

TEST_F(ChisquaredTest, InverseCovarianceMatrix) {
  ASSERT_NEAR(chisquared_.GetInverseCovarianceMatrix(0, 0), 7316.5477977290002,
              C::maxError);
}
