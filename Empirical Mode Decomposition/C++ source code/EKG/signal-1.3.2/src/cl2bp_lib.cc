// Copyright (c) 2008-2009, Evgeni A. Nurminski <nurmi@dvo.ru>
// Copyright (c) 2008-2009, Pete Gonzalez <pgonzalez@bluel.com>
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS  // disable spurious warnings
#define _USE_MATH_DEFINES // for math.h
#endif

#include "cl2bp_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#include <stdexcept>

#define BIG_NUMBER 100000

//-----------------------------------------------------------------------------------------------------------
#ifdef CL2BP_LOGGING
static void log_printf(const char *format, ...) {
  va_list argptr;
  va_start(argptr, format);

  char buf[20*1024];
  if (vsnprintf(buf,20*1024-1,format,argptr) == -1) {
    strcpy(buf, "#ERROR#");
  }
  va_end(argptr);

  cl2bp_log(buf);
}
#else
#define log_printf(...)  ((void)0)   // don't log anything (for a big speed improvement)
#endif

//-----------------------------------------------------------------------------------------------------------
static int local_max(const MallocArray<double>& x, int n, MallocArray<int>& ix) {
  int i, mx;

  mx = 0;

  MallocArray<double> xx(n+2);

  xx[0] = xx[n + 1] = -BIG_NUMBER;
  for ( i = 1; i <= n; i++)
    xx[i] = x[i - 1];
  for ( i = 0; i < n; i++ ) {
    (((xx[i] <  xx[i + 1]) && (xx[i + 1] >  xx[i + 2])) ||
     ((xx[i] <  xx[i + 1]) && (xx[i + 1] >= xx[i + 2])) ||
     ((xx[i] <= xx[i + 1]) && (xx[i + 1] >  xx[i + 2])) )
    && ( ix[mx++] = i );
  }
  return mx;
}

//-----------------------------------------------------------------------------------------------------------
static int solvps(MallocArray<double>& a, MallocArray<double>& b, int n,
  void (*cancel_handler)(void *), void *cancel_state) {

  double t;
  int j, k;
  int a_p;
  int a_q;
  int a_r;
  int a_s;

  // In empirical tests, 2% to 6% of the execution time was spent in solvps()
  if (cancel_handler) cancel_handler(cancel_state);

  for (j = 0, a_p = 0; j < n; ++j, a_p += n+1) {
    for (a_q = j*n; a_q < a_p; ++a_q)
      a[a_p] -= a[a_q] * a[a_q];
    if (a[a_p] <= 0.)
      return -1;
    a[a_p] = sqrt(a[a_p]);
    for (k = j + 1, a_q = a_p + n; k < n; ++k, a_q += n) {
      for (a_r = j*n, a_s = k*n, t = 0.; a_r < a_p;)
        t += a[a_r++] * a[a_s++];
      a[a_q] -= t;
      a[a_q] /= a[a_p];
    }
  }
  for (j = 0, a_p = 0; j < n; ++j, a_p += n+1) {
    for (k = 0, a_q = j*n; k < j;)
      b[j] -=b [k++] * a[a_q++];
    b[j] /= a[a_p];
  }
  for (j = n - 1, a_p = n*n - 1; j >= 0; --j, a_p -= n + 1) {
    for (k = j + 1, a_q = a_p + n; k < n; a_q += n)
      b[j] -= b[k++]* a[a_q];
    b[j] /= a[a_p];
  }
  return 0;
}

//-----------------------------------------------------------------------------------------------------------
static void rmmult(MallocArray<double>& rm, const MallocArray<double>& a, const MallocArray<double>& b,
  int n, int m, int l,
  void (*cancel_handler)(void *), void *cancel_state) {

  double z;
  int i,j,k;
  int b_p; // index into b
  int a_p; // index into a
  int rm_start = 0;  // index into rm
  int rm_q; // index into rm

  MallocArray<double> q0(m);
  for (i = 0; i < l ; ++i, ++rm_start) {
    // In empirical tests, 87% to 95% of the execution time was spent in rmmult()
    if (cancel_handler) cancel_handler(cancel_state);
    for (k = 0, b_p = i; k < m; b_p += l)
      q0[k++] = b[b_p];
    for (j = 0, a_p = 0, rm_q = rm_start; j < n; ++j, rm_q += l) {
      for (k = 0, z = 0.; k < m;)
        z += a[a_p++] * q0[k++];
      rm[rm_q]=z;
    }
  }
}

//-----------------------------------------------------------------------------------------------------------
static void mattr(MallocArray<double>& a, const MallocArray<double>& b, int m, int n) {
  int i, j;
  int b_p;
  int a_p = 0;
  int b_start = 0;
  for (i = 0; i < n; ++i, ++b_start)
    for ( j = 0, b_p = b_start; j < m; ++j, b_p += n )
      a[a_p++] = b[b_p];
}

//-----------------------------------------------------------------------------------------------------------
static void calcAx(MallocArray<double>& Ax, int m, int L) {
  double r = M_SQRT2, pi = M_PI;
  int i, j;

  Ax.resize((L+1)*(m + 1));

  for ( i = 0; i <= L; i++)
    for ( j = 0; j <= m; j++)
      Ax[i*(m+1) + j] = cos(i*j*pi/L);
  for ( i = 0; i < (L+1)*(m+1); i += m + 1 )
    Ax[i] /= r;
}

//-----------------------------------------------------------------------------------------------------------
#ifdef CL2BP_LOGGING
static double L2error(const MallocArray<double>& x, const MallocArray<double>& w, int L, double w1, double w2) {
  double xx, ww, sumerr, pi = M_PI;
  int i;
  for ( i = 0, sumerr = 0; i < L + 1; i++ ) {
    ww = w[i];
    xx = x[i];
    sumerr += ( ww < w1*pi || ww > w2*pi ) ?  xx*xx : (1 - xx)*(1 - xx);
  }
  return sumerr;
}
#endif // CL2BP_LOGGING
//-----------------------------------------------------------------------------------------------------------
static double CHerror(double *wmin, const MallocArray<double>& x, const MallocArray<double>& w,
  int L, double w1, double w2) {

  double xx, ww, err, errmax;
  int i;
  errmax = -1;
  *wmin = 0;
  for ( i = 0; i < L + 1; i++ ) {
    ww = w[i];
    xx = x[i];
    err = (( ww < w1 ) || (ww > w2 )) ?  fabs(xx) : fabs(1 - xx);
    if ( err > errmax ) {
      errmax = err;
      *wmin = ww;
    }
  }
  return errmax;
}

//-----------------------------------------------------------------------------------------------------------
static void Ggen(MallocArray<double>& G, int m, const MallocArray<double>& w, const MallocArray<int>& kmax,
  int l_kmax, const MallocArray<int>& kmin, int l_kmin) {

  G.resize((l_kmax + l_kmin)*(m + 1));

  int nsys, i, j;
  double r = M_SQRT2;

  nsys = l_kmax + l_kmin;
  for ( i = 0; i < l_kmax; i++)
    for ( j = 0; j <= m; j++)
      G[i*(m+1) + j] = cos(w[kmax[i]]*j);
  for ( i = l_kmax; i < nsys; i++)
    for ( j = 0; j <= m; j++)
      G[i*(m+1) + j] = -cos(w[kmin[i - l_kmax]]*j);
  for ( i = 0; i < nsys*(m+1); i += m+1 )
    G[i] /= r;
}

//-----------------------------------------------------------------------------------------------------------
bool cl2bp(MallocArray<double>& h, int m, double w1, double w2, double up[3], double lo[3], int L,
  double eps, int mit, void (*cancel_handler)(void *), void *cancel_state) {

  // Ensure sane values for input parameters
  if (m < 2 || m > 5000)
    throw std::invalid_argument("cl2bp: The m count must be between 2 and 5000");

  if (w1 < 0 || w1 > w2 || w2 > 2*M_PI)
    throw std::invalid_argument("cl2bp: The cutoffs must be in the range 0 < w1 < w2 < 2*pi");

  // Z is allocated as Z(2*L-1-2*m)
  if (L <= m)
    throw std::invalid_argument("cl2bp: The \"gridsize\" parameter must be larger than \"m\"");

  double r = M_SQRT2;
  double pi = M_PI;
  double wmin, ww, xw;
  int q1, q2, q3, i, iter = 0;
  int off, neg;

  int ik;
  int l_kmax;
  int l_kmin;
  int l_okmax;
  int l_okmin;
  double uvo = 0, lvo = 0;

  double err, diff_up, diff_lo;
  double erru, erro;
  int iup, ilo;

  int nsys, j;

  int imin = BIG_NUMBER;
  double umin;

  int k1 = -1, k2 = -1, ak1, ak2;
  double cerr;

  h.resize(2*m+1);

  bool converged = true;

  MallocArray<double> x0(L+1);
  MallocArray<double> x1(L+1);
  MallocArray<double> xx(L+1);
  MallocArray<double> xm1(L+1);
  MallocArray<double> w(L+1);
  MallocArray<double> u(L+1);
  MallocArray<double> l(L+1);
  MallocArray<double> a(L+1);
  MallocArray<double> c(L+1);
  MallocArray<int> kmax(L+1);
  MallocArray<int> kmin(L+1);
  l_kmax  = l_kmin = 0;

  MallocArray<int> okmin(L+1);
  MallocArray<int> okmax(L+1);
  l_okmax = l_okmin = 0;
  MallocArray<double> rhs(m+2);
  MallocArray<double> mu(2*(L+1));

  for ( i = 0; i <= L; i++ )
    w[i] = i*pi/L;

  MallocArray<double> Z(2*L-1-2*m);

  q1 = (int)floor(L*w1/pi);
  q2 = (int)floor(L*(w2-w1)/pi);
  q3 = L + 1 - q1 - q2;

  off = 0;
  for ( i = 0; i < q1; i++) {
    u[i] = up[0];
    l[i] = lo[0];
  }
  off += i;
  for ( i = 0; i < q2; i++) {
    u[off + i] = up[1];
    l[off + i] = lo[1];
  }
  off += i;
  for ( i = 0; i < q3; i++) {
    u[off + i] = up[2];
    l[off + i] = lo[2];
  }


  c[0] = (w2-w1)*r;
  for ( i = 1; i <= m; i++)
    c[i] =  2*(sin(w2*i)-sin(w1*i))/i;
  for ( i = 0; i <= m; i++) {
    c[i] /=  pi;
    a[i] = c[i];
  }
  log_printf("\nInitial approximation. Unconstrained L2 filter coefficients.\n");
  log_printf("=============================================================\n");

  log_printf("\nZero order term %8.3lf\n", a[0]);
  for ( i = 1; i <= m; i++) {
    log_printf("%4d %8.3lf", i, a[i]);
    if (i - 8*(i/8) == 0)
      log_printf("\n");
  }
  log_printf("\n");

  // Precalculate Ax
  MallocArray<double> Ax;
  calcAx(Ax, m, L);

  //calcA(x0, a, m, L);
  rmmult(x0, Ax, a, L + 1, m + 1, 1, cancel_handler, cancel_state);

  err = CHerror(&wmin, x0, w, L, w1, w2);
  log_printf("\nChebyshev err %12.4e (%11.5e)  <> L2 err %12.4e\n", err, wmin/pi, L2error(x0, w, L, w1, w2)/(L+1));
  for (iter = 1; ; iter++) {
    log_printf("\n:::::: Iteration %3d :::::::::::\n", iter);

    if ( (uvo > eps*2) || (lvo > eps*2) ) {
      log_printf("\nXXX Take care of old errors: uvo lvo %12.3e %12.3e", uvo, lvo);
      if( k1 >= 0 ) log_printf(" k1 %3d(%d)", k1, okmax[k1]);
      if( k2 >= 0 ) log_printf(" k2 %3d(%d)", k2, okmin[k2]);
      log_printf("\n");

      if ( uvo > lvo ) {
        kmax[l_kmax] = okmax[k1];
        l_kmax++;
        l_okmax--;
        for (i = k1; i < l_okmax; i++)
          okmax[i] = okmax[i + 1];
      } else {
        kmin[l_kmin] = okmin[k2];
        l_kmin++;
        l_okmin--;
        for (i = k2; i < l_okmin; i++)
          okmin[i] = okmin[i + 1];
      }
      nsys = l_kmax + l_kmin;

      /*
        for (i = 0; i < l_kmax; i++)
          log_printf("i %3d kmax %3d mu %12.4e\n",
             i, kmax[i], mu[i]);
        log_printf("\n");
        for (i = 0; i < l_kmin; i++)
          log_printf("i %3d kmin %3d mu %12.4e\n",
             i, kmin[i], mu[i + l_kmax]);
        log_printf("\n\n");
      */
    } else {

      //calcA(x1, a, m, L);
      rmmult(x1, Ax, a, L + 1, m + 1, 1, cancel_handler, cancel_state);


      for ( i = 0; i < l_kmax; i++)
        okmax[i] = kmax[i];
      for ( i = 0; i < l_kmin; i++)
        okmin[i] = kmin[i];
      l_okmax = l_kmax;
      l_okmin = l_kmin;

      l_kmax = local_max(x1, L + 1, kmax);


      for ( i = 0; i < L + 1; i++)
        xm1[i] = -x1[i];
      l_kmin = local_max(xm1, L + 1, kmin);

      log_printf("\nSignificant deviations from the ideal filter. Levels:");
      log_printf(" %8.2e %8.2e %8.2e (lo)", lo[0], lo[1], lo[2]);
      log_printf(" %8.2e %8.2e %8.2e (up)", up[0], up[1], up[2]);
      log_printf("\n");

      for ( i = 0, ik = 0; i < l_kmax; i++) {
        j = kmax[i];
        if ( x1[j] > u[j] - 10*eps ) {
          kmax[ik] = j;
          ik++;
        }
      }
      l_kmax = ik;

      log_printf("overshoots = %d\n", l_kmax);
      for ( i = 0; i < l_kmax; i++) {
        j = kmax[i];
        ww = w[j];
        xw = x1[j];
        err = (w1 < ww && ww < w2) ? xw - 1 : xw;
        log_printf(" i = %3d kmax = %3d xw = %12.5e err = %10.3e u = %12.4e max = %9.2e\n",
               i, j, xw, err, u[j], xw - u[j]);
      }

      for ( i = 0, ik = 0; i < l_kmin; i++) {
        j = kmin[i];
        if ( x1[j] < l[j] + 10*eps ) {
          kmin[ik] = j;
          ik++;
        }
      }
      l_kmin = ik;

      log_printf("undershoots = %d\n", l_kmin);
      for ( i = 0; i < l_kmin; i++) {
        j = kmin[i];
        ww = w[j];
        xw = -xm1[j];
        err =(w1 < ww && ww < w2) ? xw - 1 : xw;
        log_printf(" i = %3d kmin = %3d xw = %12.5e err = %10.3e l = %12.4e min = %9.2e\n",
               i, j, xw, err, l[j], xw - l[j]);
      }

      err = erru = erro = diff_up = diff_lo = 0;
      iup = ilo = 0;
      for ( i = 0; i < l_kmax; i++) {
        ik = kmax[i];
        diff_up = x1[ik] - u[ik];
        if ( diff_up > erru ) {
          erru = diff_up;
          iup = i;
        }
      }
      for ( i = 0; i < l_kmin; i++) {
        ik = kmin[i];
        diff_lo = l[ik] - x1[ik];
        if ( diff_lo > erro ) {
          erro = diff_lo;
          ilo = i;
        }
      }
      err = erro > erru ? erro : erru;
      log_printf("erru = %14.6e iup = %4d kmax = %4d ", erru, iup, kmax[iup]);
      log_printf("erro = %14.6e ilo = %4d kmin = %4d\n", erro, ilo, kmin[ilo]);
#ifndef CL2BP_LOGGING
      static_cast<void>(iup);
#endif

      if ( err < eps )
        break;
    }


    nsys = l_kmax + l_kmin;

    MallocArray<double> G(nsys*(m + 1));
    MallocArray<double> GT(nsys*(m + 1));
    MallocArray<double> GG(nsys*nsys);

    for ( i = 0; i < l_kmax; i++)
      for ( j = 0; j <= m; j++)
        G[i*(m+1) + j] = cos(w[kmax[i]]*j);
    for ( i = l_kmax; i < nsys; i++)
      for ( j = 0; j <= m; j++)
        G[i*(m+1) + j] = -cos(w[kmin[i - l_kmax]]*j);
    for ( i = 0; i < nsys*(m+1); i += m+1 )
      G[i] /= r;
    mattr(GT, G, nsys, m + 1);
    rmmult(GG, G, GT, nsys, m + 1, nsys, cancel_handler, cancel_state);


    rmmult(rhs, G, c, nsys, m + 1, 1, cancel_handler, cancel_state);
    for ( i = 0; i < nsys; i++)
      if ( i < l_kmax )
        rhs[i] -= u[kmax[i]];
      else
        rhs[i] += l[kmin[i - l_kmax]];

    for ( i = 0; i < nsys; ++i)
      mu[i] = rhs[i];

    solvps(GG, mu, nsys, cancel_handler, cancel_state);
    log_printf("\nXXX KT-system with %d equations resolved.\n", nsys);
    for ( i = 0, neg = 0; i < nsys; i++)
      if ( mu[i] < 0 ) neg++;
    log_printf("\nTotal number of negative multipliers = %3d\n\n", neg);
    while (neg) {


      umin = BIG_NUMBER;
      for ( i = 0, j = 0; i < nsys; i++) {
        if ( mu[i] >= 0 ) continue;
        j++;
        if ( mu[i] < umin ) {
          imin = i;
          umin = mu[i];
        }
      }

      if ( umin > 0 )
        break;
      log_printf(" neg = %3d of %3d imin = %3d mu-min = %12.4e\n", j, nsys, imin, umin);

      for ( i = imin; i < nsys - 1; i++) {
        rhs[i] = rhs[i + 1];
        for ( j = 0; j <= m; j++)
          G[i*(m + 1) + j] =  G[(i + 1)*(m + 1) + j];
      }

      if ( imin < l_kmax ) {
        for ( i = imin; i < l_kmax - 1; i++)
          kmax[i] = kmax[i + 1];
        l_kmax--;
      } else {
        for ( i = imin; i < nsys - 1; i++)
          kmin[i - l_kmax] = kmin[i - l_kmax + 1];
        l_kmin--;
      }

      --nsys;

      mattr(GT, G, nsys, m + 1);
      rmmult(GG, G, GT, nsys, m + 1, nsys, cancel_handler, cancel_state);
      for ( i = 0; i < nsys; ++i)
        mu[i] = rhs[i];
      solvps(GG, mu, nsys, cancel_handler, cancel_state);
    }

    MallocArray<double> da(m + 1);
    MallocArray<double> zd(nsys);

    rmmult(da, GT, mu, m + 1, nsys, 1, cancel_handler, cancel_state);
    for ( i = 0; i <= m; i++)
      a[i] = c[i] - da[i];
    rmmult(zd, G, a, nsys, m + 1, 1, cancel_handler, cancel_state);

    MallocArray<double> zz(l_okmax + l_okmin);
    Ggen(G, m, w, okmax, l_okmax, okmin, l_okmin);
    rmmult(zz, G, a, l_okmax + l_okmin, m + 1, 1, cancel_handler, cancel_state);
    uvo = lvo = eps;
    k1 = k2 = -1;
    ak1 = ak2 = -1;
    if (l_okmax + l_okmin > 0)
      log_printf("\nErrors on the previous set of freqs\n\n");
    for (i = 0; i < l_okmax; i++) {
      j = okmax[i];
      cerr = zz[i] - u[j];
      log_printf(" i %2d j %3d u %12.4e Ga %12.4e cerr %12.4e\n",
             i, j, u[j], zz[i], cerr);
      if ( cerr > uvo ) {
        uvo = cerr;
        k1 = i;
        ak1 = j;
      }
    }
    cerr = 0;
    for (i = 0; i < l_okmin; i++) {
      j = okmin[i];
      cerr = l[j] + zz[i + l_okmax];
      log_printf(" i %2d j %3d l %12.4e Ga %12.4e cerr %12.4e\n",
             i, j, l[j], zz[i + l_okmax], cerr);
      if ( cerr > lvo ) {
        lvo = cerr;
        k2 = i, ak2 = j;
      }
    }
    if ( l_okmax + l_okmin > 0 ) {
      log_printf("\n uvo = %12.4e k1 = %4d (%3d) ", uvo, k1, ak1);
      log_printf("  lvo = %12.4e k2 = %4d (%3d) ", lvo, k2, ak2);
      log_printf(" maxerr = %12.4e\n", uvo > lvo ? uvo : lvo);
#ifndef CL2BP_LOGGING
      static_cast<void>(ak1);
#endif
    }

    log_printf("\nConstrained L2 band filter coefficients.\n");
    log_printf("=====================================\n");

#ifdef CL2BP_LOGGING
    log_printf("\nZero order term %8.3lf\n", a[0]);
    for ( i = 1; i <= m; i++) {
      log_printf("%4d %8.3lf", i, a[i]);
      if (i - 8*(i/8) == 0)
        log_printf("\n");
    }
    log_printf("\n");
#endif

    //calcA(xx, a, m, L);
    rmmult(xx, Ax, a, L + 1, m + 1, 1, cancel_handler, cancel_state);

    if ( iter >= mit ) {
      log_printf("Maximum iterations reached\n");
      converged = false;
    }
  }
  for (i = 0; i < m; i++) {
    h[i] = a[m - i]/2;
    h[m + i + 1] = a[i + 1]/2;
  }
  h[m] = a[0]*r/2;

  return converged;
}
