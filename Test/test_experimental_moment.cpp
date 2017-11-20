// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../Src/experimental_moment.h"
#include "../Src/Weights.h"

const double maxError = 1e-15;

typedef ExperimentalMoment E;
typedef Weights W;

TEST(WRatioTest, Matthias) {
  E experiment;
  ASSERT_NEAR(experiment.wRatio(3., W::WD00, 1), 1, maxError);
}
