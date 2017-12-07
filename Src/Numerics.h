// Copyright 2017 Dirk Hornung

#ifndef PROGRAM_SRC_NUMERICS_H_
#define PROGRAM_SRC_NUMERICS_H_

#include <vector>
#include <cmath>
#include "nr3.h"
#include "./constants.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef Constants C;

namespace Ublas {

  using boost::numeric::ublas::matrix;
  using boost::numeric::ublas::permutation_matrix;
  using boost::numeric::ublas::identity_matrix;
  using std::vector;

// Numerical Operations
// -> Integration : Gaussian Quadratures
class Numerics {
 public:
  template <class Iter>
  double kahan_summation(Iter begin, Iter end) {
    double result = 0.;

    double c = 0.;
    for (; begin != end; ++begin) {
      double y = *begin - c;
      double t = result + y;
      c = (t - result) - y;
      result = t;
    }
    return result;
  }

  /* Matrix inversion routine.
  Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
  template<class T>
  static bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse) {
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<T> A(input);

    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    // perform LU-factorization
    int res = lu_factorize(A, pm);
    if (res != 0)
      return false;

    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<T> (A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    return true;
  }

  // Returns weights and abcisas needed for Gaussian Quadratures
  template <typename T>
    static void gauleg(const T x1, const T x2, std::vector<T> &x,
                       std::vector<T> &w) {
    const double EPS = 1.0e-15;
    T z1, z, xm, xl, pp, p3, p2, p1;
    int n = x.size();
    int m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);

    for (int i=0; i < m; i++) {
      z = cos(3.141592654*(i+0.75)/(n+0.5));

      do {
        p1 = 1.0;
        p2 = 0.0;
        for (int j=0; j < n; j++) {
          p3 = p2;
          p2 = p1;
          T jComplx(j);
          p1 = ((2.0*jComplx+1.0)*z*p2-jComplx*p3)/(jComplx+1.);
        }
        T nComplx(n);
        pp = nComplx*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1-p1/pp;
      } while (abs(z-z1) > EPS);
      x[i] = xm-xl*z;
      x[n-1-i] = xm+xl*z;
      w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
      w[n-1-i] = w[i];
    }
  }

 public:
  template <class T>
  static void outputVector(vector<T> vector) {
    std::cout << "Output Vector:" << std::endl;
    for (int i = 0; i < 2; i++) {
      std::cout << vector[i] << "\t";
    }
    std::cout << std::endl << std::endl;
  }

  template <class T>
  static void outputMatrix(matrix<T> matrix) {
    // int rows = matrix.size1()
    // int cols = matrix.size2()
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        std::cout << matrix(i, j) << "\t";
      }
      std::cout << std::endl;
    }
    std::cout << matrix.size1() << std::endl;
  }
  
  // Integration: Gaussian Quadratures
  // Integrates any function from a to b
  template <typename T, typename Func>
  static T Integrate(Func &func, const T a, const T b) {
    vector<T> w(1201);
    vector<T> x(1201);
    Numerics::gauleg(a, b, x, w);
    T s = 0;
    for(int j=0; j<x.size(); j++) {
      s += w[j]*func(x[j]);
    }
    return s;
  };

  static void gaussj(MatDoub_IO &a, MatDoub_IO &b) {
    Int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
    Doub big,dum,pivinv;
    VecInt indxc(n),indxr(n),ipiv(n);
    for (j=0;j<n;j++) ipiv[j]=0;
    for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++)
        if (ipiv[j] != 1)
          for (k=0;k<n;k++) {
            if (ipiv[k] == 0) {
              if (abs(a[j][k]) >= big) {
                big=abs(a[j][k]);
                irow=j;
                icol=k;
              }
            }
          }
      ++(ipiv[icol]);
      if (irow != icol) {
        for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
        for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=0;l<n;l++) a[icol][l] *= pivinv;
      for (l=0;l<m;l++) b[icol][l] *= pivinv;
      for (ll=0;ll<n;ll++)
        if (ll != icol) {
          dum=a[ll][icol];
          a[ll][icol]=0.0;
          for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
          for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
        }
    }
    for (l=n-1;l>=0;l--) {
      if (indxr[l] != indxc[l])
        for (k=0;k<n;k++)
          SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
  }

  static void gaussj(MatDoub_IO &a) {
    MatDoub b(a.nrows(), 0);
    gaussj(a,b);
  }
};

}  // namespace Ublas

using Ublas::Numerics;

#endif
