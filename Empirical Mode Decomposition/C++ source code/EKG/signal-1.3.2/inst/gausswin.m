## Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
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
## @deftypefn  {Function File} {} gausswin (@var{m})
## @deftypefnx {Function File} {} gausswin (@var{m}, @var{a})
##
## Return the filter coefficients of a Gaussian window of length @var{m}.
## The width of the window is inversely proportional to the parameter @var{a}.
## Use larger @var{a} for a narrow window.  Use larger @var{m} for a smoother
## curve.
##
##     w = exp ( -(a*x)^2/2 )
##
## for x = linspace(-(m-1)/m, (m-1)/m, m)
## @end deftypefn

function w = gausswin (m, a)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("gausswin: M must be a positive integer");
  elseif (nargin == 1)
    a = 2.5;
  endif

  w = exp ( -0.5 * ( a/m * [ -(m-1) : 2 : m-1 ]' ) .^ 2 );

endfunction

%!assert (gausswin (1), 1)

%% Test input validation
%!error gausswin ()
%!error gausswin (0.5)
%!error gausswin (-1)
%!error gausswin (ones (1, 4))
%!error gausswin (1, 2, 3)
