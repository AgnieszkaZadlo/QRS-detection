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
## @deftypefn  {Function File} {@var{n} =} cheb2ord (@var{wp}, @var{ws}, @var{rp}, @var{rs})
## @deftypefnx {Function File} {@var{n} =} cheb2ord ([@var{wp1}, @var{wp2}], [@var{ws1}, @var{ws2}], @var{rp}, @var{rs})
## @deftypefnx {Function File} {[@var{n}, @var{wc}] =} cheb2ord (@dots{})
## Compute the minimum filter order of a Chebyshev type II filter with the
## desired response characteristics. The filter frequency band edges are
## specified by the passband frequency @var{wp} and stopband frequency @var{ws}.
## Frequencies are normalized to the Nyquist frequency in the range [0,1].
## @var{rp} is the allowable passband ripple measured in decibels, and @var{rs}
## is the minimum attenuation in the stop band, also in decibels.  The output
## arguments @var{n} and @var{wc} can be given as inputs to @code{cheby2}.
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
## @seealso{buttord, cheb1ord, cheby2, ellipord}
## @end deftypefn

function [n, Wc] = cheb2ord(Wp, Ws, Rp, Rs)

  if nargin != 4
    print_usage;
  else
    validate_filter_bands ("cheb2ord", Wp, Ws);
  endif

  T = 2;

  ## returned frequency is the same as the input frequency
  Wc = Ws;

  ## warp the target frequencies according to the bilinear transform
  Ws = (2/T)*tan(pi*Ws./T);
  Wp = (2/T)*tan(pi*Wp./T);

  if (Wp(1) < Ws(1))
    ## low pass
    if (length(Wp) == 1)
      Wa = Wp/Ws;
    else
      ## FIXME: Implement band reject filter type
      error ("cheb2ord: band reject is not yet implemented");
    endif;
  else
   ## if high pass, reverse the sense of the test
   if (length(Wp) == 1)
      Wa = Ws/Wp;
    else
      ## band pass
      Wa=(Wp.^2 - Ws(1)*Ws(2))./(Wp*(Ws(1)-Ws(2)));
    endif;
  endif;
  Wa = min(abs(Wa));

  ## compute minimum n which satisfies all band edge conditions
  stop_atten = 10^(abs(Rs)/10);
  pass_atten = 10^(abs(Rp)/10);
  n = ceil(acosh(sqrt((stop_atten-1)/(pass_atten-1)))/acosh(1/Wa));

endfunction

%!demo
%! Fs = 10000;
%! [n, Wc] = cheb2ord (1000/(Fs/2), 1200/(Fs/2), 0.5, 29);
%!
%! subplot (221);
%! plot ([0, 1000, 1000, 0, 0], [0, 0, -0.5, -0.5, 0], ";;");
%! hold on;
%! grid;
%! title("Pass band Wp=1000 Rp=0.0");
%! xlabel("Frequency (Hz)");
%! ylabel("Attenuation (dB)");
%! [b, a] = cheby2 (n, 29, Wc);
%! [h, w] = freqz (b, a, [], Fs);
%! plot (w, 20*log10(abs(h)), ";;");
%! axis ([ 0, 1500, -1, 0]);
%! hold off;
%!
%! subplot (222);
%! plot ([1200, Fs/2, Fs/2, 1200, 1200], [-29, -29, -500, -500, -29], ";;");
%! hold on;
%! axis ([ 0, Fs/2, -250, 0]);
%! title("Stop band Ws=1200 Rs=29");
%! xlabel("Frequency (Hz)");
%! ylabel("Attenuation (dB)");
%! grid;
%! [b, a] = cheby2 (n, 29, Wc);
%! [h, w] = freqz (b, a, [], Fs);
%! plot (w, 20*log10(abs(h)), ";;");
%! hold off;
%!
%! subplot (223);
%! plot ([0, 1000, 1000, 0, 0], [0, 0, -0.5, -0.5, 0], ";;");
%! hold on;
%! axis ([ 800, 1010, -0.6, -0.0]);
%! title("Pass band detail Wp=1000 Rp=0.5");
%! xlabel("Frequency (Hz)");
%! ylabel("Attenuation (dB)");
%! grid;
%! [b, a] = cheby2 (n, 29, Wc);
%! [h, w] = freqz (b, a, [800:1010], Fs);
%! plot (w, 20*log10(abs(h)), "r;filter n;");
%! [b, a] = cheby2 (n-1, 29, Wc);
%! [h, w] = freqz (b, a, [800:1010], Fs);
%! plot (w, 20*log10(abs(h)), "b;filter n-1;");
%! [b, a] = cheby2 (n+1, 29, Wc);
%! [h, w] = freqz (b, a, [800:1010], Fs);
%! plot (w, 20*log10(abs(h)), "g;filter n+1;");
%! hold off;
%!
%! subplot (224);
%! plot ([1200, Fs/2, Fs/2, 1200, 1200], [-29, -29, -500, -500, -29], ";;");
%! hold on;
%! axis ([ 1190, 1210, -40, -20]);
%! title("Stop band detail Wp=1200 Rp=29");
%! xlabel("Frequency (Hz)");
%! ylabel("Attenuation (dB)");
%! grid;
%! [b, a] = cheby2 (n, 29, Wc);
%! [h, w] = freqz (b, a, [1190:1210], Fs);
%! plot (w, 20*log10(abs(h)), "r;filter n;");
%! [b, a] = cheby2 (n-1, 29, Wc);
%! [h, w] = freqz (b, a, [1190:1210], Fs);
%! plot (w, 20*log10(abs(h)), "b;filter n-1;");
%! [b, a] = cheby2 (n+1, 29, Wc);
%! [h, w] = freqz (b, a, [1190:1210], Fs);
%! plot (w, 20*log10(abs(h)), "g;filter n+1;");
%! hold off;

%% Test case from demo
%!assert (cheb2ord (0.2, 0.24, 0.5, 29), 8)

%% Test input validation
%!error cheb2ord ()
%!error cheb2ord (.1)
%!error cheb2ord (.1, .2)
%!error cheb2ord (.1, .2, 3)
%!error cheb2ord (.1, .2, 3, 4, 5)
%!error cheb2ord ([.1 .1], [.2 .2], 3, 4)
%!error cheb2ord ([.1 .2], [.5 .6], 3, 4)
%!error cheb2ord ([.1 .5], [.2 .6], 3, 4)
