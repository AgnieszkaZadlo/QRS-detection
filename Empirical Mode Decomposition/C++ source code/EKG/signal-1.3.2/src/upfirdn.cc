// Copyright (C) 2008 Eric Chassande-Mottin, CNRS (France) <ecm@apc.univ-paris7.fr>
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

#include <octave/config.h>
#include <octave/defun-dld.h>
#include <octave/error.h>
#include <octave/gripes.h>
#include <octave/oct-obj.h>
#include <octave/pager.h>
#include <octave/quit.h>
#include <octave/variables.h>

template<class MT>
MT upfirdn (MT &x, ColumnVector &h, octave_idx_type p, octave_idx_type q)
{
  octave_idx_type rx = x.rows ();
  octave_idx_type cx = x.columns ();
  bool isrowvector = false;

  if ((rx == 1) && (cx > 1)) // if row vector, transpose to column vector
    {
      x = x.transpose ();
      rx = x.rows ();
      cx = x.columns ();
      isrowvector = true;
    }

  octave_idx_type Lh = h.length ();
  const double r = p/(static_cast<double> (q));

  const octave_idx_type Ly = ceil (static_cast<double> ((rx-1)*p + Lh) /
                                   static_cast<double> (q));

  MT y (Ly, cx, 0.0);

  for (octave_idx_type c = 0; c < cx; c++)
    {
      octave_idx_type m = 0;
      while (m < Ly)
        {
          const octave_idx_type n = floor (m/r);
          const octave_idx_type lm = (m * q) % p;
          octave_idx_type k = 0;
          typename MT::element_type accum;
          accum = 0.0;
          do
            {
              octave_idx_type ix = n - k;
              if (ix >= rx)
                {
                  k ++;
                  continue;
                }

              const octave_idx_type ih = k * p + lm;
              if ((ih >= Lh) | (ix < 0))
                break;

              accum += h (ih) * x (ix, c);
              k++;
            }
          while (1);

          y (m, c) = accum;
          m ++;
        }

    }

  if (isrowvector)
    y = y.transpose ();

  return y;
}

DEFUN_DLD (upfirdn, args,,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{y} =} upfirdn (@var{x}, @var{h}, @var{p}, @var{q})\n\
Upsample, FIR filtering, and downsample.\n\
@end deftypefn\n")
{
  octave_value_list retval;

  const int nargin = args.length ();

  if (nargin != 4)
    {
      print_usage ();
      return retval;
    }

  ColumnVector h (args (1).vector_value ());

  if (error_state)
    {
      gripe_wrong_type_arg ("upfirdn", args (1));
      return retval;
    }

  octave_idx_type p = args (2).idx_type_value ();

  if (error_state)
    {
      gripe_wrong_type_arg ("upfirdn", args (2));
      return retval;
    }

  octave_idx_type q = args (3).idx_type_value ();

  if (error_state)
    {
      gripe_wrong_type_arg ("upfirdn", args (3));
      return retval;
    }

  // Do the dispatching
  if (args (0).is_real_type ())
    {
      Matrix x = args (0).matrix_value ();
      if (error_state)
        {
          gripe_wrong_type_arg ("upfirdn", args (0));
          return retval;
        }

      Matrix y = upfirdn (x, h, p, q);
      retval (0) = y;
    }
  else if (args (0).is_complex_type ())
    {
      ComplexMatrix x = args (0).complex_matrix_value ();
      if (error_state)
        {
          gripe_wrong_type_arg ("upfirdn", args (0));
          return retval;
        }

      ComplexMatrix y = upfirdn (x, h, p, q);
      retval (0) = y;
    }
  else
    {
      gripe_wrong_type_arg ("upfirdn", args (0));
      return retval;
    }

  return retval;
}

/*
%!assert (isequal (upfirdn (1:100, 1, 1, 1), 1:100))
%!assert (isequal (upfirdn (1:100, 1, 1, 2), 1:2:100))

%% Test input validation
%!error upfirdn ()
%!error upfirdn (1,2)
%!error upfirdn (1,2,3)
%!error upfirdn (1,2,3,4,5)
*/

