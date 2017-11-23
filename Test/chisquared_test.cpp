// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../Src/chisquared.h"

const double maxError = 1e-13;

TEST(InvertedCovarianceMatrix, Matthias) {
  ASSERT_NEAR(1., 1., maxError);
}
