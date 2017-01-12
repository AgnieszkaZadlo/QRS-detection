## Copyright (C) 2014 Mike Miller <mtmiller@ieee.org>
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
## @deftypefn  {Function File} {} validate_filter_bands (@var{func}, @var{wp}, @var{ws})
## Validate filter design frequency response bands.  This is an internal
## convenience function used by @code{buttord}, @code{cheb1ord},
## @code{cheb2ord}, @code{ellipord}.
## @end deftypefn

function validate_filter_bands (func, wp, ws)

  if (nargin != 3)
    print_usage ();
  elseif (! ischar (func))
    error ("validate_filter_bands: FUNC must be a string");
  elseif (! (isvector (wp) && isvector (ws) && (numel (wp) == numel (ws))))
    error ([func ": WP and WS must both be scalars or vectors of length 2\n"]);
  elseif (! ((numel (wp) == 1) || (numel (wp) == 2)))
    error ([func ": WP and WS must both be scalars or vectors of length 2\n"]);
  elseif (! (isreal (wp) && all (wp >= 0) && all (wp <= 1)))
    error ([func ": all elements of WP must be in the range [0,1]\n"]);
  elseif (! (isreal (ws) && all (ws >= 0) && all (ws <= 1)))
    error ([func ": all elements of WS must be in the range [0,1]\n"]);
  elseif ((numel (wp) == 2) && (wp(2) <= wp(1)))
    error ([func ": WP(1) must be less than WP(2)\n"])
  elseif ((numel (ws) == 2) && (ws(2) <= ws(1)))
    error ([func ": WS(1) must be less than WS(2)\n"])
  elseif ((numel (wp) == 2) && (all (wp > ws) || all (ws > wp)))
    error ([func ": WP must be contained by WS or WS must be contained by WP\n"]);
  endif

endfunction
