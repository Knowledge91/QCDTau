#include <iostream>
#include "gtest/gtest.h"

TEST(IndependentMethod, ResetsToZero) {
  int i = 3;
  EXPECT_EQ(0, i);
}
