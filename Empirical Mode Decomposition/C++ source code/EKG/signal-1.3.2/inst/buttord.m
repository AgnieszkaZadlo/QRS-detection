## Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
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
## @deftypefn  {Function File} {@var{n} =} buttord (@var{wp}, @var{ws}, @var{rp}, @var{rs})
## @deftypefnx {Function File} {@var{n} =} buttord ([@var{wp1}, @var{wp2}], [@var{ws1}, @var{ws2}], @var{rp}, @var{rs})
## @deftypefnx {Function File} {[@var{n}, @var{wc}] =} buttord (@dots{})
## Compute the minimum filter order of a Butterworth filter with the desired
## response characteristics.  The filter frequency band edges are specified by
## the passband frequency @var{wp} and stopband frequency @var{ws}.  Frequencies
## are normalized to the Nyquist frequency in the range [0,1].  @var{rp} is the
## allowable passband ripple measured in decibels, and @var{rs} is the minimum
## attenuation in the stop band, also in decibels.  The output arguments @var{n}
## and @var{wc} can be given as inputs to @code{butter}.
##
## If @var{wp} and @var{ws} are scalars, then @var{wp} is the passband cutoff
## frequency and @var{ws} is the stopband edge frequency.  If @var{ws} is
## greater than @var{wp}, the filter is a low-pass filter.  If @var{wp} is
## greater than @var{ws}, the filter is a high-pass filter.
##
## If @var{wp} and @var{ws} are vectors of length 2, then @var{wp} defines the
## passband interval and @var{ws} defines the stopband interval.  If @var{wp}
## is contained within @var{ws} (@var{ws1} < @var{wp1} < @var{wp2} < @var{ws2}),
## the filter is a band-pass filter.  If @var{ws} is contained within @var{wp}
## (@var{wp1} < @var{ws1} < @var{ws2} < @var{wp2}), the filter is a band-stop
## or band-reject filter.
##
## Theory: |H(W)|^2 = 1/[1+(W/Wc)^(2N)] = 10^(-R/10)
## With some algebra, you can solve simultaneously for Wc and N given
## Ws,Rs and Wp,Rp.  For high pass filters, subtracting the band edges
## from Fs/2, performing the test, and swapping the resulting Wc back
## works beautifully.  For bandpass and bandstop filters this process
## significantly overdesigns.  Artificially dividing N by 2 in this case
## helps a lot, but it still overdesigns.
##
## @seealso{butter, cheb1ord, cheb2ord, ellipord}
## @end deftypefn

function [n, Wc] = buttord(Wp, Ws, Rp, Rs)

  if (nargin != 4)
    print_usage ();
  else
    validate_filter_bands ("buttord", Wp, Ws);
  endif

  if (numel (Wp) == 2)
    warning ("buttord: seems to overdesign bandpass and bandreject filters");
  endif

  T = 2;

  ## if high pass, reverse the sense of the test
  stop = find(Wp > Ws);
  Wp(stop) = 1-Wp(stop); # stop will be at most length 1, so no need to
  Ws(stop) = 1-Ws(stop); # subtract from ones(1,length(stop))

  ## warp the target frequencies according to the bilinear transform
  Ws = (2/T)*tan(pi*Ws./T);
  Wp = (2/T)*tan(pi*Wp./T);

  ## compute minimum n which satisfies all band edge conditions
  ## the factor 1/length(Wp) is an artificial correction for the
  ## band pass/stop case, which otherwise significantly overdesigns.
  qs = log(10^(Rs/10) - 1);
  qp = log(10^(Rp/10) - 1);
  n = ceil(max(0.5*(qs - qp)./log(Ws./Wp))/length(Wp));

  ## compute -3dB cutoff given Wp, Rp and n
  Wc = exp(log(Wp) - qp/2/n);

  ## unwarp the returned frequency
  Wc = atan(T/2*Wc)*T/pi;

  ## if high pass, reverse the sense of the test
  Wc(stop) = 1-Wc(stop);

endfunction

%% Test input validation
%!error buttord ()
%!error buttord (.1)
%!error buttord (.1, .2)
%!error buttord (.1, .2, 3)
%!error buttord (.1, .2, 3, 4, 5)
%!error buttord ([.1 .1], [.2 .2], 3, 4)
%!error buttord ([.1 .2], [.5 .6], 3, 4)
%!error buttord ([.1 .5], [.2 .6], 3, 4)
