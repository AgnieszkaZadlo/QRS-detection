## Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
## Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
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
## @deftypefn  {Function File} {[@var{b}, @var{a}] =} cheby1 (@var{n}, @var{rp}, @var{w})
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} cheby1 (@var{n}, @var{rp}, @var{w}, "high")
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} cheby1 (@var{n}, @var{rp}, [@var{wl}, @var{wh}])
## @deftypefnx {Function File} {[@var{b}, @var{a}] =} cheby1 (@var{n}, @var{rp}, [@var{wl}, @var{wh}], "stop")
## @deftypefnx {Function File} {[@var{z}, @var{p}, @var{g}] =} cheby1 (@dots{})
## @deftypefnx {Function File} {[@var{a}, @var{b}, @var{c}, @var{d}] =} cheby1 (@dots{})
## @deftypefnx {Function File} {[@dots{}] =} cheby1 (@dots{}, "s")
## Generate a Chebyshev type I filter with @var{rp} dB of passband ripple.
##
## [b, a] = cheby1(n, Rp, Wc)
##    low pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby1(n, Rp, Wc, 'high')
##    high pass filter with cutoff pi*Wc radians
##
## [b, a] = cheby1(n, Rp, [Wl, Wh])
##    band pass filter with edges pi*Wl and pi*Wh radians
##
## [b, a] = cheby1(n, Rp, [Wl, Wh], 'stop')
##    band reject filter with edges pi*Wl and pi*Wh radians
##
## [z, p, g] = cheby1(...)
##    return filter as zero-pole-gain rather than coefficients of the
##    numerator and denominator polynomials.
##
## [...] = cheby1(...,'s')
##     return a Laplace space filter, W can be larger than 1.
##
## [a,b,c,d] = cheby1(...)
##  return  state-space matrices
##
## References:
##
## Parks & Burrus (1987). Digital Filter Design. New York:
## John Wiley & Sons, Inc.
## @end deftypefn

function [a, b, c, d] = cheby1 (n, rp, w, varargin)

  if (nargin > 5 || nargin < 3 || nargout > 4 || nargout < 2)
    print_usage ();
  endif

  ## interpret the input parameters
  if (! (isscalar (n) && (n == fix (n)) && (n > 0)))
    error ("cheby1: filter order N must be a positive integer");
  endif

  stop = false;
  digital = true;
  for i = 1:numel (varargin)
    switch (varargin{i})
      case "s"
        digital = false;
      case "z"
        digital = true;
      case {"high", "stop"}
        stop = true;
      case {"low", "pass"}
        stop = false;
      otherwise
        error ("cheby1: expected [high|stop] or [s|z]");
    endswitch
  endfor

  if (! ((numel (w) <= 2) && (rows (w) == 1 || columns (w) == 1)))
    error ("cheby1: frequency must be given as WC or [WL, WH]");
  elseif ((numel (w) == 2) && (w(2) <= w(1)))
    error ("cheby1: W(1) must be less than W(2)");
  endif

  if (digital && ! all ((w >= 0) & (w <= 1)))
    error ("cheby1: all elements of W must be in the range [0,1]");
  elseif (! digital && ! all (w >= 0))
    error ("cheby1: all elements of W must be in the range [0,inf]");
  endif

  if (! (isscalar (rp) && isnumeric (rp) && (rp >= 0)))
    error ("cheby1: passband ripple RP must be a non-negative scalar");
  endif

  ## Prewarp to the band edges to s plane
  if (digital)
    T = 2;       # sampling frequency of 2 Hz
    w = 2 / T * tan (pi * w / T);
  endif

  ## Generate splane poles and zeros for the Chebyshev type 1 filter
  C = 1;  ## default cutoff frequency
  epsilon = sqrt (10^(rp / 10) - 1);
  v0 = asinh (1 / epsilon) / n;
  pole = exp (1i * pi * [-(n - 1):2:(n - 1)] / (2 * n));
  pole = -sinh (v0) * real (pole) + 1i * cosh (v0) * imag (pole);
  zero = [];

  ## compensate for amplitude at s=0
  gain = prod (-pole);
  ## if n is even, the ripple starts low, but if n is odd the ripple
  ## starts high. We must adjust the s=0 amplitude to compensate.
  if (rem (n, 2) == 0)
    gain = gain / 10^(rp / 20);
  endif

  ## splane frequency transform
  [zero, pole, gain] = sftrans (zero, pole, gain, w, stop);

  ## Use bilinear transform to convert poles to the z plane
  if (digital)
    [zero, pole, gain] = bilinear (zero, pole, gain, T);
  endif

  ## convert to the correct output form
  if (nargout == 2)
    a = real (gain * poly (zero));
    b = real (poly (pole));
  elseif (nargout == 3)
    a = zero;
    b = pole;
    c = gain;
  else
    ## output ss results
    [a, b, c, d] = zp2ss (zero, pole, gain);
  endif

endfunction

%% Test input validation
%!error [a, b] = cheby1 ()
%!error [a, b] = cheby1 (1)
%!error [a, b] = cheby1 (1, 2)
%!error [a, b] = cheby1 (1, 2, 3, 4, 5, 6)
%!error [a, b] = cheby1 (.5, 2, .2)
%!error [a, b] = cheby1 (3, 2, .2, "invalid")

