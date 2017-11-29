// Copyright 2017 Dirk Hornung

#include <vector>
#include "gtest/gtest.h"
#include "../src/constants.h"
#include "../src/chisquared.h"


const double maxError = 1e-13;

//Constants constants(3, 3, 4);
std::vector<double> s0_set = {C::sTau, 3., 2.8, 2.6, 2.4, 2.3, 2.2, 2.1, 2.};
//Chisquared chisquared(constants, s0_set);

TEST(InvertedCovarianceMatrix, Matthias) {
  //  ASSERT_NEAR(chisquared.GetInvertedCovarianceMatrix(0, 0),
  //          10848.945884521445, maxError);
}
