// Copyright 2017 Dirk Hornung

#include <functional>
#include "./experimental_moment.h"
#include "./experimental_data.h"
#include "./weights.h"

using std::function;

ExperimentalMoment::ExperimentalMoment(
    double s0, function<double(double)> weight) :
    s0_(s0), weight_(weight), data_(ExperimentalData()) {}

int ExperimentalMoment::getBinNumber() {
    double pos = 0;
    for (int i=0; i < 80; i++) {
      pos += data_.GetDSbin(i);
      if (pos >= s0_) {
        return i+1;
      }
    }
    return 80;
}

