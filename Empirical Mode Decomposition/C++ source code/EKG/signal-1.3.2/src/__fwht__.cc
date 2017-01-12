// Copyright (C) 2013 Mike Miller <mtmiller@ieee.org>
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

#include <octave/oct.h>

template <typename T>
static void do_fwht (T *x, octave_idx_type n)
{
  double *xa = &x[0];
  double *xb = &x[n/2];

  // Perform the FWHT butterfly on the input
  for (octave_idx_type i = 0; i < n/2; i++)
    {
      double ya = xa[i] + xb[i];
      double yb = xa[i] - xb[i];
      xa[i] = ya;
      xb[i] = yb;
    }

  // Divide and conquer
  if (n > 2)
    {
      do_fwht (xa, n/2);
      do_fwht (xb, n/2);
    }
}

DEFUN_DLD (__fwht__, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} __fwht__ (@var{x})\n\
Undocumented internal function.\n\
@end deftypefn")
{
  int nargin = args.length ();
  if (nargin != 1)
    print_usage ();
  else
    {
      NDArray x = args(0).array_value ();
      octave_idx_type n = x.rows ();
      double *data = x.fortran_vec ();
      for (octave_idx_type i = 0; i < x.columns (); i++)
        {
          do_fwht (&data[i*n], n);
        }
      return octave_value (x);
    }
  return octave_value ();
}

/*
## No test needed for internal helper function.
%!assert (1)
*/
