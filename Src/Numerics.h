// Numerics.h
#ifndef NUMERICS_H
#define NUMERICS_H

#include <cmath>
#include <vector>

using namespace std;


// Numerical Operations
// -> Integration : Gaussian Quadratures
class Numerics {
  // Returns weights and abcisas needed for Gaussian Quadratures (Numerics::qgauss)
  template <typename T>
  static void gauleg(const T x1, const T x2, vector<T> &x, vector<T> &w) {
    const double EPS=1.0e-15;
    T z1, z, xm, xl, pp, p3, p2, p1;
    int n = x.size();
    int m = (n+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);

    for (int i=0; i<m; i++) {
      z = cos( 3.141592654*(i+0.75)/(n+0.5));

      do {
        p1 = 1.0;
        p2 = 0.0;
        for (int j=0; j<n; j++) {
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
      w[n-1-i]=w[i];
    }
  }

 public:
  // Integration: Gaussian Quadratures
  // Integrates any function from a to b
  template <typename T, typename Func>
  static T qgauss(Func &func, const T a, const T b) {
    vector<T> w(1201);
    vector<T> x(1201);
    Numerics::gauleg(a, b, x, w);
    T s = 0;
    for(int j=0; j<x.size(); j++) {
      s += w[j]*func(x[j]);
    }
    return s;
  };

};
#endif
