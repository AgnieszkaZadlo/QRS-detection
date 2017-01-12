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

#ifndef CL2BP_H
#define CL2BP_H

#include <cstdlib>
#include <cassert>
#include <cstring> //for memset

//-----------------------------------------------------------------------------------------------------------
// If you want to debug the cl2bp algorithm, define the CL2BP_LOGGING symbol and provide an
// implementation of cl2bp_log().
#ifdef CL2BP_LOGGING
extern void cl2bp_log(const char *message);
#endif

//-----------------------------------------------------------------------------------------------------------
// This is a simple helper class that performs bounds-checking for debug builds, but reduces to an unchecked
// malloc() array for release builds.
template <class T>
class MallocArray {
  int length;
  T *ptr;
public:
  T& operator[](int index) {
    assert(index >= 0 && index < length);
    return ptr[index];
  }
  T operator[](int index) const {
    assert(index >= 0 && index < length);
    return ptr[index];
  }

  int get_length() const { return length; }

  void resize(int length_) {
    assert(length_ >= 0 && length_ <= 512*1024*1024);  // verify that the array size is reasonable
    length = length_;
    ptr = (T *)realloc(ptr, length * sizeof(T));
    std::memset(ptr, 0, length * sizeof(T));
  }

  MallocArray(int length_=0) {
    ptr = 0;
    length = 0;
    if (length_) resize(length_);
  }
  ~MallocArray() { free(ptr); }
private:
  MallocArray(const MallocArray&) { }  // copy constructor is unimplemented and disallowed
};

//-----------------------------------------------------------------------------------------------------------
// Constrained L2 BandPass FIR filter design
//  h       2*m+1 filter coefficients (output)
//  m       degree of cosine polynomial
//  w1,w2   fist and second band edges
//  up      upper bounds
//  lo      lower bounds
//  L       grid size
//  eps     stopping condition ("SN")
//  mit     maximum number of iterations
//
// cl2bp() returns true if the solver converged, or false if the maximum number of iterations was reached.
// If provided, the cancel_handler function pointer will be called periodically inside long-running loops,
// giving an opportunity to abort (by throwing a C++ exception).  The cancel_state parameter is a
// user-defined pointer, e.g. a button object, boolean flag, or other means of detecting a cancel request.
//
// Example usage:
//   MallocArray<double> coefficients;
//   double up[3] = { 0.02, 1.02, 0.02 };
//   double lo[3] = { -0.02, 0.98, -0.02 };
//   cl2bp(coefficients, 30, 0.3*M_PI, 0.6*M_PI, up, lo, 1<<11, 1.e-5, 20);
//
// The algorithm is described in this paper:
//   I. W. Selesnick, M. Lang, and C. S. Burrus.  A modified algorithm for constrained least square
//   design of multiband FIR filters without specified transition bands.  IEEE Trans. on Signal
//   Processing, 46(2):497-501, February 1998.
bool cl2bp(MallocArray<double>& h, int m, double w1, double w2, double up[3], double lo[3], int L,
  double eps, int mit, void (*cancel_handler)(void *)=0, void *cancel_state=0);

#endif
