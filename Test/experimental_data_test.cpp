// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../src/experimental_data.h"


const double maxError = 1e-13;

ExperimentalData data;

TEST(AlephData, Matthias) {
  ASSERT_NEAR(data.GetSbin(0), 0.03749999999999999, maxError);
  ASSERT_NEAR(data.GetDErr(0), 4.6889399699999999e-4, maxError);
}

TEST(ErrorMatrix, Matthias) {
  ASSERT_NEAR(data.GetErrorMatrix(0, 0), 2.1986158042263601e-7, maxError);
  ASSERT_NEAR(data.GetErrorMatrix(80, 80), 1.6000000000000001e-3, maxError);
  ASSERT_NEAR(data.GetErrorMatrix(81, 81), 3.7134137767198075e-5, maxError);
}
