## Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
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
## @deftypefn {Function File} {} boxcar (@var{m})
## Return the filter coefficients of a rectangular window of length @var{m}.
## @end deftypefn

function w = boxcar (m)

  if (nargin != 1)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("boxcar: M must be a positive integer");
  endif

  w = ones(m, 1);

endfunction

%!assert (boxcar (1), 1)
%!assert (boxcar (2), ones (2, 1))
%!assert (boxcar (100), ones (100, 1))

%% Test input validation
%!error boxcar ()
%!error boxcar (0.5)
%!error boxcar (-1)
%!error boxcar (ones (1, 4))
%!error boxcar (1, 2)
