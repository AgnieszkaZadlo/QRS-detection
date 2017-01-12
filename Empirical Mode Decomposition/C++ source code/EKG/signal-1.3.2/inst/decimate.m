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
## @deftypefn  {Function File} {@var{y} =} decimate (@var{x}, @var{q})
## @deftypefnx {Function File} {@var{y} =} decimate (@var{x}, @var{q}, @var{n})
## @deftypefnx {Function File} {@var{y} =} decimate (@dots{}, "fir")
##
## Downsample the signal @var{x} by a reduction factor of @var{q}. A lowpass
## antialiasing filter is applied to the signal prior to reducing the input
## sequence. By default, an order @var{n} Chebyshev type I filter is used.
## If @var{n} is not specified, the default is 8.
##
## If the optional argument @code{"fir"} is given, an order @var{n} FIR filter
## is used, with a default order of 30 if @var{n} is not given.
##
## Note that @var{q} must be an integer for this rate change method.
##
## Example:
## @example
## ## Generate a signal that starts away from zero, is slowly varying
## ## at the start and quickly varying at the end, decimate and plot.
## ## Since it starts away from zero, you will see the boundary
## ## effects of the antialiasing filter clearly.  Next you will see
## ## how it follows the curve nicely in the slowly varying early
## ## part of the signal, but averages the curve in the quickly
## ## varying late part of the signal.
## t = 0:0.01:2;
## x = chirp (t, 2, .5, 10, "quadratic") + sin (2*pi*t*0.4);
## y = decimate (x, 4);
## stem (t(1:121) * 1000, x(1:121), "-g;Original;"); hold on; # original
## stem (t(1:4:121) * 1000, y(1:31), "-r;Decimated;"); hold off; # decimated
## @end example
## @end deftypefn

function y = decimate(x, q, n, ftype)

  if (nargin < 2 || nargin > 4)
    print_usage ();
  elseif (! (isscalar (q) && (q == fix (q)) && (q > 0)))
    error ("decimate: Q must be a positive integer");
  endif

  if (nargin < 3)
    ftype = "iir";
    n = [];
  elseif (nargin < 4)
    if (ischar (n))
      ftype = n;
      n = [];
    else
      ftype = "iir";
    endif
  endif

  if (! any (strcmp (ftype, {"fir", "iir"})))
    error ('decimate: filter type must be either "fir" or "iir"');
  endif

  fir = strcmp (ftype, "fir");
  if (isempty (n))
    if (fir)
      n = 30;
    else
      n = 8;
    endif
  endif

  if (! (isscalar (n) && (n == fix (n)) && (n > 0)))
    error ("decimate: N must be a positive integer");
  endif

  if (fir)
    b = fir1 (n, 1/q);
    y = fftfilt (b, x);
  else
    [b, a] = cheby1 (n, 0.05, 0.8/q);
    y = filtfilt (b, a, x);
  endif

  y = y(1:q:length(x));

endfunction

%!demo
%! t=0:0.01:2; x=chirp(t,2,.5,10,'quadratic')+sin(2*pi*t*0.4);
%! y = decimate(x,4);   # factor of 4 decimation
%! stem(t(1:121)*1000,x(1:121),"-g;Original;"); hold on; # plot original
%! stem(t(1:4:121)*1000,y(1:31),"-r;Decimated;"); hold off; # decimated
%! %------------------------------------------------------------------
%! % The signal to decimate starts away from zero, is slowly varying
%! % at the start and quickly varying at the end, decimate and plot.
%! % Since it starts away from zero, you will see the boundary
%! % effects of the antialiasing filter clearly.  You will also see
%! % how it follows the curve nicely in the slowly varying early
%! % part of the signal, but averages the curve in the quickly
%! % varying late part of the signal.

%% Test input validation
%!error decimate ()
%!error decimate (1)
%!error decimate (1, 2, 3, 4, 5)
%!error decimate (1, -1)

