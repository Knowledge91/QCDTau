// Copyright 2017 Dirk Hornung

#include "gtest/gtest.h"
#include "../Src/experimental_moment.h"
#include "../Src/Weights.h"

const double maxError = 1e-13;

typedef ExperimentalMoment E;
typedef Weights W;

E experiment(3.1570893124000001, W::WD00);


TEST(AlephData, Matthias) {
  ASSERT_NEAR(experiment.sbin[0], 0.03749999999999999, maxError);
}

TEST(WRatioTest, Matthias) {
  ASSERT_NEAR(experiment.wRatio(2.1, W::WD00, 0),
              0.99930452590445074, maxError);
}

TEST(Jacobian, Matthias) {
  ASSERT_NEAR(experiment.jacobian[0],
              0.056094687833062200, maxError);
}
