## Copyright (C) 2013 - Juan Pablo Carbajal
##
## This progrm is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{y} =} upsamplefill (@var{x}, @var{v})
## @deftypefnx {Function File} {@var{y} =} upsamplefill (@dots{}, @var{copy})
## Upsamples a vector interleaving given values or copies of the vector elements.
##
## The values in the vector @var{v} are placed between the elements of @var{x}.
##
## If the optional argument @var{copy} is @var{true} then @var{v} should be a
## scalar and each value in @var{x} are repeat @var{v} times.
##
## Example:
## @example
## @group
## upsamplefill (eye (2), 2, true)
## @result{}  1   0
##     1   0
##     1   0
##     0   1
##     0   1
##     0   1
## upsamplefill (eye (2), [-1 -1 -1])
## @result{}  1   0
##    -1  -1
##    -1  -1
##    -1  -1
##     0   1
##    -1  -1
##    -1  -1
##    -1  -1
## @end group
## @end example
##
## @seealso{upsample}
## @end deftypefn

## Author: Juan Pablo Carbajal <ajuanpi+dev@gmail.com>

function y = upsamplefill (x, v, copy=false)

  if nargin<2
    print_usage;
  end

  [nr,nc] = size (x);
  if copy

    if any ([nr,nc]==1)

      y = kron (x(:), ones(v+1,1));
      if nr == 1
        y = y.';
      endif

    else

      y = kron (x, ones(v+1,1));

    endif

    return

  else

    % Assumes 'v' row or column vector
    n = length(v) + 1;
    N = n*nr;

    if any ([nr,nc]==1)

      N        = N*nc;
      idx      = 1:n:N;
      idx_c    = setdiff (1:N, 1:n:N);
      y        = zeros (N,1);
      y(idx)   = x;
      y(idx_c) = repmat (v(:), max(nr,nc), 1);

      if nr == 1
        y = y.';
      endif

    else

      idx      = 1:n:N;
      idx_c    = setdiff(1:N,1:n:N);
      y        = zeros (N,nc);
      y(idx,:)   = x;

      y(idx_c,:) = repmat (v(:), nr, nc);

    endif
  endif

endfunction

%!assert(upsamplefill([1,3,5],2),[1,2,3,2,5,2]);
%!assert(upsamplefill([1;3;5],2),[1;2;3;2;5;2]);
%!assert(upsamplefill([1,2,5],[2 -2]),[1,2,-2,2,2,-2,5,2,-2]);
%!assert(upsamplefill(eye(2),2,true),[1,0;1,0;1,0;0,1;0,1;0,1]);
%!assert(upsamplefill([1,3,5],2,true),[1,1,1,3,3,3,5,5,5]);
%!assert(upsamplefill([1;3;5],2,true),[1;1;1;3;3;3;;5;5;5]);
