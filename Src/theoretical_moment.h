// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_THEORETICAL_MOMENT_H_
#define PROGRAM_SRC_THEORETICAL_MOMENT_H_

#include <cmath>
#include <complex>
#include "./constants.h"
#include "./numerics.h"
#include "CRunDec.h"


class TheoreticalMoment {
 public:
  explicit TheoreticalMoment(Constants constants) : constants_(constants) {
    c_run_dec_ = new CRunDec(constants.getNf());
  }

  std::complex<double> GetSpectralMoment(double s0);

  std::complex<double> GetD0(std::complex<double> s);

  std::complex<double> GetL(std::complex<double> s, double mu);

  double GetAlphaMu(double mu);

  Constants constants_;
  CRunDec* c_run_dec_;
};

#endif  // PROGRAM_SRC_THEORETICAL_MOMENT_H_
