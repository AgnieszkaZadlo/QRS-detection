## Author: Paul Kienzle <pkienzle@users.sf.net> (2004)
## This program is granted to the public domain.

## -*- texinfo -*-
## @deftypefn  {Function File} {} flattopwin (@var{m})
## @deftypefnx {Function File} {} flattopwin (@var{m}, "periodic")
## @deftypefnx {Function File} {} flattopwin (@var{m}, "symmetric")
##
## Return the filter coefficients of a Flat Top window of length @var{m}.
## The Flat Top window is defined by the function f(w):
##
## @example
## @group
##   f(w) = 1 - 1.93 cos(2 pi w) + 1.29 cos(4 pi w)
##            - 0.388 cos(6 pi w) + 0.0322cos(8 pi w)
## @end group
## @end example
##
## where w = i/(m-1) for i=0:m-1 for a symmetric window, or
## w = i/m for i=0:m-1 for a periodic window.  The default
## is symmetric.  The returned window is normalized to a peak
## of 1 at w = 0.5.
##
## This window has low pass-band ripple, but high bandwidth.
##
## According to [1]:
##
##    The main use for the Flat Top window is for calibration, due
##    to its negligible amplitude errors.
##
## [1] Gade, S; Herlufsen, H; (1987) "Use of weighting functions in DFT/FFT
## analysis (Part I)", Bruel & Kjaer Technical Review No.3.
## @end deftypefn

function w = flattopwin (m, opt)

  if (nargin < 1 || nargin > 2)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("flattopwin: M must be a positive integer");
  endif

  N = m - 1;
  if (nargin == 2)
    switch (opt)
      case "periodic"
        N = m;
      case "symmetric"
        N = m - 1;
      otherwise
        error ("flattopwin: window type must be either \"periodic\" or \"symmetric\"");
    endswitch
  endif

  if (m == 1)
    w = 1;
  else
    x = 2*pi*[0:m-1]'/N;
    w = (1-1.93*cos(x)+1.29*cos(2*x)-0.388*cos(3*x)+0.0322*cos(4*x))/4.6402;
  endif

endfunction

%!assert (flattopwin (1), 1);
%!assert (flattopwin (2), 0.0042 / 4.6402 * ones (2, 1), eps);
%!assert (flattopwin (15), flipud (flattopwin (15)), 10*eps);
%!assert (flattopwin (16), flipud (flattopwin (16)), 10*eps);

%!assert (flattopwin (15), flattopwin (15, "symmetric"));
%!assert (flattopwin (16)(1:15), flattopwin (15, "periodic"));

%% Test input validation
%!error flattopwin ()
%!error flattopwin (0.5)
%!error flattopwin (-1)
%!error flattopwin (ones (1, 4))
%!error flattopwin (1, 2)
%!error flattopwin (1, 2, 3)
%!error flattopwin (1, "invalid")

