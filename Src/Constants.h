// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_CONSTANTS_H_
#define PROGRAM_SRC_CONSTANTS_H_

#include <cmath>
#include <vector>

class Constants {
 public:
  Constants(int nc, int nf, int loops);

  int getNf() { return nf_; }
  int getNc() { return nc_; }
  double getZeta(int i) { return zeta_[i]; }
  double getBeta(int i) { return beta_[i]; }
  double getC(int i, int j) { return c_[i][j]; }

  const double kAlphaSMz = 0.1181;  // PDG 2016 p.145
  const double kMz = 91.1876;  // PDG 2016 p.29
  static constexpr double kMu = 1.77682;  // sqrt(sTau)
  static constexpr double kAlphaSTau = 0.31927;
  static constexpr double kSTau = 3.1570893124;  // mTau*mTau;

  static constexpr double sTau = 3.1570893124;  // mTau*mTau;
  static constexpr double Be = 17.827;         // HFAG 2011
  static constexpr double kDBe = 0.04;          // HFAG 2011
  static constexpr double kRtauVex = 1.;
  static constexpr double kDRtauVex = 0.;
  static constexpr double mu2 = 1.;
  static constexpr double mTau = 1.77686;       // PDG 2016

  static constexpr double Vud = 0.97425;        // Towner, Hardy 2009
  static constexpr double dVud = 0.00022;
  static constexpr double SEW = 1.0198;         // EW radiative corr.
  static constexpr double dSEW = 0.0006;

  // Pseudoscalar resonance parameter
  static constexpr double fpi = 92.21e-3;       // PDG 2010
  static constexpr double dfpi = 0.14e-3;       // PDG 2010

  // particle masses
  static constexpr double mpim = 0.13957018;    // M_pi^-

  // Pi-Pole
  static double pifac;
  static double dpifac;

  static constexpr double maxError = 1e-15;

 private:
  double loops_;  // number of loops
  double nc_;  // colour number
  double nf_;  // flavour number
  double beta_[5];  // beta coefficients
  double zeta_[8];  // zeta function
  double c_[6][6];  // Adler coefficients
};

#endif  // PROGRAM_SRC_CONSTANTS_H_
