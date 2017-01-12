## Copyright (C) 2014 Mike Miller <mtmiller@ieee.org>
##
## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} hann (@var{m})
## @deftypefnx {Function File} {} hann (@var{m}, "periodic")
## @deftypefnx {Function File} {} hann (@var{m}, "symmetric")
## Return the filter coefficients of a Hanning window of length @var{m}.
##
## If the optional argument @code{"periodic"} is given, the periodic form
## of the window is returned.  This is equivalent to the window of length
## @var{m}+1 with the last coefficient removed.  The optional argument
## @code{"symmetric"} is equivalent to not specifying a second argument.
##
## This function exists for @sc{matlab} compatibility only, and is equivalent
## to @code{hanning (@var{m})}.
##
## @seealso{hanning}
## @end deftypefn

function w = hann (varargin)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  w = hanning (varargin{:});

endfunction

%!assert (hann (1), 1);
%!assert (hann (2), zeros (2, 1));
%!assert (hann (16), flipud (hann (16)), 10*eps);
%!assert (hann (15), flipud (hann (15)), 10*eps);
%!test
%! N = 15;
%! A = hann (N);
%! assert (A(ceil (N/2)), 1);

%!assert (hann (15), hann (15, "symmetric"));
%!assert (hann (16)(1:15), hann (15, "periodic"));
%!test
%! N = 16;
%! A = hann (N, "periodic");
%! assert (A (N/2 + 1), 1);

%% Test input validation
%!error hann ()
%!error hann (0.5)
%!error hann (-1)
%!error hann (1, "invalid")

