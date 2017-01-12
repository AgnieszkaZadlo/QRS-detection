/*
 * Copyright 2000 Paul Kienzle <pkienzle@users.sf.net>
 * This source code is freely redistributable and may be used for
 * any purpose.  This copyright notice must be maintained.
 * Paul Kienzle is not responsible for the consequences of using
 * this software.
 *
 * Mar 2000 - Kai Habel (kahacjde@linux.zrz.tu-berlin.de)
 *      Change: ColumnVector x=arg(i).vector_value();
 *      to: ColumnVector x=ColumnVector(arg(i).vector_value());
 * Oct 2000 - Paul Kienzle (pkienzle@users.sf.net)
 *      rewrite to ignore NaNs rather than replacing them with zero
 *      extend to handle matrix arguments
 */

#include <cmath>
#include <octave/oct.h>
#include <octave/lo-ieee.h>
#include <octave/pager.h>
using namespace std;

// The median class holds a sorted data window.  This window is
// intended to slide over the data, so when the window shifts
// by one position, the old value from the start of the window must
// be removed and the new value from the end of the window must
// be added.  Since removals and additions generally occur in pairs,
// a hole is left in the sorted window when the value is removed so
// that on average, fewer values need to be shifted to close the
// hole and open a new one in the sorted position.
class Median {
private:
  double *window; // window data
  int max;        // length of window used
  int hole;       // position of hole, or max if no hole
  void close_hole() // close existing hole
  {
    // move hole to the end of the window
    while (hole < max-1) {
      window[hole] = window[hole+1];
      hole++;
    }
    // shorten window (if no hole, then hole==max)
    if (hole == max-1) max--;
  }
  void print();

public:
  Median(int n)
  {
    max=hole=0;
    window = new double[n];
  }
  ~Median(void)
  {
    delete [] window;
  }

  void add(double v);          // add a new value
  void remove(double v);       // remove an existing value
  void clear() { max=hole=0; } // clear the window
  double operator() ();        // find the median in the window
} ;

// Print the sorted window, and indicate any hole
void Median::print()
{
  octave_stdout << "[ ";
  for (int i=0; i < max; i++)
    {
      if (i == hole)
        octave_stdout << "x ";
      else
        octave_stdout << window[i] << " ";
    }
  octave_stdout << " ]";
}


// Remove a value from the sorted window, leaving a hole.  The caller
// must promise to only remove values that they have added.
void Median::remove(double v)
{
  // NaN's are not added or removed
  if (lo_ieee_isnan(v)) return;

  //  octave_stdout << "Remove " << v << " from "; print();

  // only one hole allowed, so close pre-existing ones
  close_hole();

  // binary search to find the value to remove
  int lo = 0, hi=max-1;
  hole = hi/2;
  while (lo <= hi) {
    if (v > window[hole]) lo = hole+1;
    else if (v < window[hole]) hi = hole-1;
    else break;
    hole = (lo+hi)/2;
  }

  // Verify that it is the correct value to replace
  // Note that we shouldn't need this code since we are always replacing
  // a value that is already in the window, but for some reason
  // v==window[hole] occasionally doesn't work.
  if (v != window[hole]) {
    for (hole = 0; hole < max-1; hole++)
      if (fabs(v-window[hole]) < fabs(v-window[hole+1])) break;
    warning ("medfilt1: value %f not found---removing %f instead",
             v, window[hole]);
    print(); octave_stdout << endl;
  }

  //  octave_stdout << " gives "; print(); octave_stdout << endl;
}

// Insert a new value in the sorted window, plugging any holes, or
// extending the window as necessary.  The caller must promise not
// to add more values than median was created with, without
// removing some beforehand.
void Median::add(double v)
{
  // NaN's are not added or removed
  if (lo_ieee_isnan(v)) return;

  //  octave_stdout << "Add " << v << " to "; print();

  // If no holes, extend the array
  if (hole == max) max++;

  // shift the hole up to the beginning as far as it can go.
  while (hole > 0 && window[hole-1] > v) {
    window[hole] = window[hole-1];
    hole--;
  }

  // or shift the hole down to the end as far as it can go.
  while (hole < max-1 && window[hole+1] < v) {
    window[hole] = window[hole+1];
    hole++;
  }

  // plug in the replacement value
  window[hole] = v;

  // close the hole
  hole = max;

  //  octave_stdout << " gives "; print(); octave_stdout << endl;
}

// Compute the median value from the sorted window
// Return the central value if there is one or the average of the
// two central values.  Return NaN if there are no values.
double Median::operator()()
{
  close_hole();

  if (max % 2 == 1)
    return window[(max-1)/2];
  else if (max == 0)
    return lo_ieee_nan_value();
  else
    return (window[max/2-1]+window[max/2])/2.0;
}

DEFUN_DLD (medfilt1, args, ,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{y} =} medfilt1 (@var{x})\n\
@deftypefnx {Loadable Function} {@var{y} =} medfilt1 (@var{x}, @var{n})\n\
Apply a median filter of length @var{n} to the signal @var{x}.  A sliding\n\
window is applied to the data, and for each step the median value in the\n\
window is returned.  If @var{n} is odd then the window for y(i) is\n\
x(i-(n-1)/2:i+(n-1)/2).  If @var{n} is even then the window is\n\
x(i-n/2:i+n/2-1) and the two values in the center of the sorted window are\n\
averaged.  If @var{n} is not given, then 3 is used.  NaNs are ignored,\n\
as are values beyond the ends, by taking the median of the remaining\n\
values.\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length();
  if (nargin < 1 || nargin > 3)
    {
      print_usage ();
      return retval;
    }

  if (args(0).is_complex_type())
    {
      error("medfilt1 cannot process complex vectors");
      return retval;
    }

  int n=3;    // length of the filter (default 3)
  if (nargin > 1) n = NINT(args(1).double_value());
  if (n < 1)
    {
      error ("medfilt1 filter length must be at least 1");
      return retval;
    }

  // Create a window to hold the sorted median values
  Median median(n);
  int mid = n/2;             // mid-point of the window

  Matrix signal(args(0).matrix_value());
  int nr = signal.rows();    // number of points to process
  int nc = signal.columns(); // number of points to process
  Matrix filter(nr,nc);      // filtered signal to return

  if (nr == 1) // row vector
    {
      int start = -n, end = 0, pos=-(n-mid)+1;
      while (pos < nc)
        {
          if (start >= 0) median.remove(signal(0,start));
          if (end < nc)   median.add(signal(0,end));
          if (pos >= 0)   filter(0,pos) = median();
          start++, end++, pos++;
        }
    }
  else // column vector or matrix
    {
      for (int column=0; column < nc; column++)
        {
          median.clear();
          int start = -n, end = 0, pos=-(n-mid)+1;
          while (pos < nr)
            {
              if (start >= 0) median.remove(signal(start,column));
              if (end < nr)   median.add(signal(end,column));
              if (pos >= 0)   filter(pos,column) = median();
              start++, end++, pos++;
            }
        }
    }

  retval(0) = filter;
  return retval;
}

/*
%!assert(medfilt1([1, 2, 3, 4, 5], 3), [1.5, 2, 3, 4, 4.5]);
*/
