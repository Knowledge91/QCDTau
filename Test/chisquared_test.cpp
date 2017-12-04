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




TEST(ChisquaredTest, CovarianceMatrix) {
  //  ASSERT_NEAR(chisquared.,
  //          10848.945884521445, maxError);
}
