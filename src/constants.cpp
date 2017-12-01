// Copyright 2017 Dirk Hornung

#include "./constants.h"


Constants::Constants(int nc, int nf, int loops) :
    nc_(nc), nf_(nf), loops_(loops) {
  zeta_[3] = 1.20205690315959;  // coefficients.nb
  zeta_[5] = 1.03692775514337;  // coefficients.nb
  zeta_[7] = 1.00834927738192;  // coefficients.nb

  beta_[1] = 11./2. - 1./3.*nf_;  // rgm06
  beta_[2] = 51./4. - 19./12.*nf_;  // rgm06
  beta_[3] = 2857./64. - 5033./576.*nf_ + 325./1728.*pow(nf_, 2);  // rgm06
  beta_[4] = 149753./768. + 891./32.*zeta_[3]  // rgm06
      -(1078361./20736. + 1627./864.*zeta_[3])*nf_
      + (50065./20736. + 809./1296.*zeta_[3])*pow(nf_, 2)
      + 1093./93312.*pow(nf_, 3);

  c_[0][0] = -5./3.; c_[0][1] = 1;  // rgm06
  c_[1][1] = 1.; c_[1][2] = 0.;  //  rgm06
  c_[2][1] = 365./24. - 11.*zeta_[3] - (11./12. - 2./3.*zeta_[3])*nf_;
  c_[2][2] = -beta_[1]*c_[1][1]/4.; c_[2][3] = 0.;  // rgm06
  c_[3][1] = 87029./288. - 1103./4.*zeta_[3] + 275./6.*zeta_[5]
      +(-7847./216. + 262./9.*zeta_[3] - 25./9.*zeta_[5])*nf_
      + (151./162.-19./27.*zeta_[3])*pow(nf_, 2);  // rgm06
  c_[3][2] = -1./4.*(beta_[2]*c_[1][1]+2*beta_[1]*c_[2][1]);
  c_[3][3] = pow(beta_[1], 2)/12.*c_[1][1]; c_[3][4] = 0.;  // rgm06
  c_[4][1] = 78631453./20736. - 1704247./432.*zeta_[3]
      + 4185./8.*pow(zeta_[3], 2) + 34165./96.*zeta_[5]
      - 1995./16.*zeta_[7];  // Diogo PHD
  c_[4][2] = -1./4.*(beta_[3]*c_[1][1]+2*beta_[2]*c_[2][1]+3*beta_[1]*c_[3][1]);
  c_[4][3] = beta_[1]/24.*(5.*beta_[2]*c_[1][1]+6*beta_[1]*c_[2][1]);  // rgm-6
  c_[4][4] = -pow(beta_[1], 3)/32.*c_[1][1]; c_[4][5] = 0.;  // rgm06
  c_[5][1] = 283.;
  c_[5][2] = 1./4.*(-beta_[4]*c_[1][1] - 2.*beta_[3]*c_[2][1]-3.
                    *beta_[2]*c_[3][1]-4.*beta_[1]*c_[4][1]);
  c_[5][3] = 1./24.*(12.*c_[3][1]*pow(beta_[1], 2)+6.*beta_[1]*beta_[3]*c_[1][1]
                  +14.*beta_[2]*beta_[1]*c_[2][1]+3.*pow(beta_[2], 2)*c_[1][1]);
  c_[5][4] = 1./96.*(-12*pow(beta_[1], 3)*c_[2][1]
                    -13.*beta_[2]*pow(beta_[1], 2)*c_[1][1]);
  c_[5][5] = 1./80.*pow(beta_[1], 4)*c_[1][1];
}

double Constants::pifac = 24.*std::pow(M_PI*Vud*fpi, 2.)*SEW;
double Constants::dpifac = Constants::pifac*std::sqrt(4.*std::pow(dVud/Vud, 2)
                         + std::pow(dSEW/SEW, 2)+4.*std::pow(dfpi/fpi, 2));




