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
## @deftypefn  {Function File} {@var{y} =} digitrevorder (@var{x}, @var{r})
## @deftypefnx {Function File} {[@var{y}, @var{i}] =} digitrevorder (@var{x}, @var{r})
## Reorder the elements of the vector @var{x} in digit-reversed order.
## The elements of @var{x} are converted to radix @var{r} and reversed.
## The reordered indices of the elements of @var{x} are returned in @var{i}.
## @seealso{bitrevorder, fft, ifft}
## @end deftypefn

function [y, i] = digitrevorder (x, r)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! isvector (x))
    error ("digitrevorder: X must be a vector");
  elseif (!(isscalar (r) && r == fix (r) && r >= 2 && r <= 36))
    error ("digitrevorder: R must be an integer between 2 and 36");
  else
    tmp = log (numel (x)) / log (r);
    if (fix (tmp) != tmp)
      error ("digitrevorder: X must have length equal to an integer power of %d", r);
    endif
  endif

  old_ind = 0:numel (x) - 1;
  new_ind = base2dec (fliplr (dec2base (old_ind, r)), r);

  i = new_ind + 1;
  y(old_ind + 1) = x(i);

  if (columns (x) == 1)
    y = y';
  else
    i = i';
  endif

endfunction

%!assert (digitrevorder (0, 2), 0);
%!assert (digitrevorder (0, 36), 0);
%!assert (digitrevorder (0:3, 4), 0:3);
%!assert (digitrevorder ([0:3]', 4), [0:3]');
%!assert (digitrevorder (0:7, 2), [0 4 2 6 1 5 3 7]);
%!assert (digitrevorder (0:15, 2), [0 8 4 12 2 10 6 14 1 9 5 13 3 11 7 15]);
%!assert (digitrevorder (0:15, 4), [0 4 8 12 1 5 9 13 2 6 10 14 3 7 11 15]);

%% Test input validation
%!error digitrevorder ();
%!error digitrevorder (1);
%!error digitrevorder (1, 2, 3);
%!error digitrevorder ([], 1);
%!error digitrevorder ([], 37);
%!error digitrevorder (0:3, 8);
