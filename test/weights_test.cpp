// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../src/weights.h"
#include "../src/constants.h"

typedef Weights W;
typedef Constants C;

TEST(WeightsTest, WD00) {
  ASSERT_NEAR(W::WD00(2.), -3. , C::maxError);
}

TEST(WeightsTest, WTau) {
  //  ASSERT_NEAR(W::WTau(2.75.), -3. , C::maxError);
}
