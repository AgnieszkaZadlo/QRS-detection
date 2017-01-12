## Copyright (C) 2007 R.G.H. Eschauzier <reschauzier@yahoo.com>
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
## @deftypefn {Function File} {@var{p_out} =} polyrev (@var{p_in})
## Undocumented internal function.  This function is used by the impinvar
## and invimpinvar functions in the signal package.
## @end deftypefn

## Adapted by Carnë Draug on 2011 <carandraug+dev@gmail.com>

## Reverse the coefficients of a polynomial

function p_out = polyrev (p_in)

  p_out = p_in(end:-1:1);

endfunction
