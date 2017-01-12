## Copyright (C) 2013   Lukas F. Reichlin
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {@var{y} =} wconv (@var{type}, @var{x}, @var{f})
## @deftypefnx {Function File} {@var{y} =} wconv (@var{type}, @var{x}, @var{f}, @var{shape})
## 1-D or 2-D convolution.
## 
## @strong{Inputs}
## @table @var
## @item type
## Type of convolution.
## @item x
## Signal vector or matrix.
## @item f
## Coefficients of @acronym{FIR} filter.
## @item shape
## Shape.
## @end table
##
## @strong{Outputs}
## @table @var
## @item y
## Convoluted signal.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2013
## Version: 0.1

function y = wconv (type, x, f, shape = "full")

  if (nargin < 3 || nargin > 4)
    print_usage ();
  endif

  switch (type(1))
    case {1, "1"}
      y = conv2 (x(:).', f(:).', shape);
      if (rows (x) > 1)
        y = y.';
      endif
    case {2, "2"}
      y = conv2 (x, f, shape);
    case "r"
      y = conv2 (x, f(:).', shape);
    case "c"
      y = conv2 (x.', f(:).', shape);
      y = y.';
    otherwise
      print_usage ();
  endswitch

endfunction
