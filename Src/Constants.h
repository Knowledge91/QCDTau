// Constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>

class Constants {
 public:
 Constants(int nf) : nf(nf) {
    zeta[3] = 1.20205690315959; // coefficients.nb
    zeta[5] = 1.03692775514337; // coefficients.nb

    beta[1] = 11./2. - 1./3.*nf; // rgm06
    beta[2] = 51./4. - 19./12.*nf; // rgm06
    beta[3] = 2857./64. - 5033./576.*nf + 325./1728.*pow(nf, 2); // rgm06
    beta[4] = 149753./768. + 891./32.*zeta[3] - (1078361./20736. + 1627./864.*zeta[3])*nf // rgm06
            + (50065./20736. + 809./1296.*zeta[3])*pow(nf,2) + 1093./93312.*pow(nf,3);

    c[0][0] = -5./3.; c[0][1] = 1; // rgm06
    c[1][1] = 1.; c[1][2] = 0.; // rgm06
    c[2][1] = 365./24. - 11.*zeta[3] - (11./12. - 2./3.*zeta[3])*nf; c[2][2] = beta[1]*c[1][1]/4.; c[2][3] = 0.; // rgm06
    c[3][1] = 87029./288. - 1103./4.*zeta[3] + 275./6.*zeta[5] + (-7847./216. + 262./9.*zeta[3] - 25./9.*zeta[5])*nf + (151./162.-19./27.*zeta[3])*pow(nf,2); // rgm06
    c[3][2] = -1./4.*(beta[2]*c[1][1]+2*beta[1]*c[2][1]); c[3][3] = pow(beta[1],2)/12.*c[1][1]; c[3][4] = 0.;
    c[4][2] = -1./4.*(beta[3]*c[1][1] + 2*beta[2]*c[2][1] + 3*beta[1]*c[3][1]); c[4][3] = beta[1]/24.*(5.*beta[2]*c[1][1] + 6*beta[1]*c[2][1]);
    c[4][4] = -pow(beta[1],3)/32.*c[1][1]; c[4][5] = 0.;

  }
  double nf; // flavour number
  double beta[5]; // beta coefficients
  double zeta[6]; // zeta function
  double c[5][6]; // Adler coefficients
  static constexpr double Be = 17.827; // HFAG 2011
  static constexpr double mu2 = 1.;
  static constexpr double mTau = 1.77686; // PDG 2016
  static constexpr double sTau = 3.1570893124; //mTau*mTau;

  static constexpr double maxError = 1e-15; // maximal allowed errro
};

#endif
