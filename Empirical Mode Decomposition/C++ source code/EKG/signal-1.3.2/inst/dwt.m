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
## @deftypefn  {Function File} {[@var{u}, @var{v}] =} dwt (@var{x}, @var{wname})
## @deftypefnx {Function File} {[@var{u}, @var{v}] =} dwt (@var{x}, @var{Hp}, @var{Gp})
## @deftypefnx {Function File} {[@var{u}, @var{v}] =} dwt (@var{x}, @var{Hp}, @var{Gp}, @dots{})
## Discrete wavelet transform (1D).
## 
## @strong{Inputs}
## @table @var
## @item x
## Signal vector.
## @item wname
## Wavelet name.
## @item Hp
## Coefficients of low-pass decomposition @acronym{FIR} filter.
## @item Gp
## Coefficients of high-pass decomposition @acronym{FIR} filter.
## @end table
##
## @strong{Outputs}
## @table @var
## @item u
## Signal vector of average, approximation.
## @item v
## Signal vector of difference, detail.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: April 2013
## Version: 0.1

function [u, v] = dwt (x, varargin)

  if (nargin == 2)
    wname = varargin{1};
    [Hp, Gp] = wfilters (wname, "d");
  elseif (nargin == 3)
    Hp = varargin{1};
    Gp = varargin{2};
  else
    print_usage ();
  endif

  tmp = wconv (1, x, Hp, "valid");
  u = tmp(1:2:end);
  
  tmp = wconv (1, x, Gp, "valid");
  v = tmp(1:2:end);

endfunction
