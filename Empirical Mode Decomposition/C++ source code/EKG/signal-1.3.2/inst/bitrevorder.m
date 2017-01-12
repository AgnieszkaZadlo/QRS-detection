## Copyright (C) 2007 Sylvain Pelissier <sylvain.pelissier@gmail.com>
## Copyright (C) 2013 Mike Miller
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
## @deftypefn  {Function File} {@var{y} =} bitrevorder (@var{x})
## @deftypefnx {Function File} {[@var{y} @var{i}] =} bitrevorder (@var{x})
## Reorder the elements of the vector @var{x} in bit-reversed order.
## Equivalent to calling @code{digitrevorder (@var{x}, 2)}.
## @seealso{digitrevorder, fft, ifft}
## @end deftypefn

function [y, i] = bitrevorder (x)

  if (nargin != 1)
    print_usage ();
  elseif (! isvector (x))
    error ("bitrevorder: X must be a vector");
  elseif (fix (log2 (numel (x))) != log2 (numel (x)))
    error ("bitrevorder: X must have length equal to an integer power of 2");
  endif

  [y, i] = digitrevorder (x, 2);

endfunction

%!assert (bitrevorder (0), 0);
%!assert (bitrevorder (0:1), 0:1);
%!assert (bitrevorder ([0:1]'), [0:1]');
%!assert (bitrevorder (0:7), [0 4 2 6 1 5 3 7]);
%!assert (bitrevorder (0:15), [0 8 4 12 2 10 6 14 1 9 5 13 3 11 7 15]);

%% Test input validation
%!error bitrevorder ();
%!error bitrevorder (1, 2);
%!error bitrevorder ([]);
%!error bitrevorder (0:2);
