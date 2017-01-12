## Copyright (C) 2001 Paulo Neis <p_neis@yahoo.com.br>
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
## @deftypefn  {Function File} {@var{n} =} ellipord (@var{wp}, @var{ws}, @var{rp}, @var{rs})
## @deftypefnx {Function File} {@var{n} =} ellipord ([@var{wp1}, @var{wp2}], [@var{ws1}, @var{ws2}], @var{rp}, @var{rs})
## @deftypefnx {Function File} {[@var{n}, @var{wc}] =} ellipord (@dots{})
## Compute the minimum filter order of an elliptic filter with the desired
## response characteristics.  The filter frequency band edges are specified
## by the passband frequency @var{wp} and stopband frequency @var{ws}.
## Frequencies are normalized to the Nyquist frequency in the range [0,1].
## @var{rp} is the allowable passband ripple measured in decibels, and @var{rs}
## is the minimum attenuation in the stop band, also in decibels.  The output
## arguments @var{n} and @var{wc} can be given as inputs to @code{ellip}.
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
## Reference: Lamar, Marcus Vinicius, @cite{Notas de aula da disciplina TE 456 -
## Circuitos Analogicos II}, UFPR, 2001/2002.
## @seealso{buttord, cheb1ord, cheb2ord, ellip}
## @end deftypefn

function [n, Wp] = ellipord(Wp, Ws, Rp, Rs)

  if (nargin != 4)
    print_usage ();
  else
    validate_filter_bands ("ellipord", Wp, Ws);
  endif

  ## sampling frequency of 2 Hz
  T = 2;

  Wpw = tan(pi.*Wp./T); # prewarp
  Wsw = tan(pi.*Ws./T); # prewarp

  ## pass/stop band to low pass filter transform:
  if (length(Wpw)==2 && length(Wsw)==2)
    wp=1;
    w02 = Wpw(1) * Wpw(2);      # Central frequency of stop/pass band (square)
    w3 = w02/Wsw(2);
    w4 = w02/Wsw(1);
    if (w3 > Wsw(1))
      ws = (Wsw(2)-w3)/(Wpw(2)-Wpw(1));
    elseif (w4 < Wsw(2))
      ws = (w4-Wsw(1))/(Wpw(2)-Wpw(1));
    else
      ws = (Wsw(2)-Wsw(1))/(Wpw(2)-Wpw(1));
    endif
  elseif (Wpw > Wsw)
    wp = Wsw;
    ws = Wpw;
  else
    wp = Wpw;
    ws = Wsw;
  endif

  k=wp/ws;
  k1=sqrt(1-k^2);
  q0=(1/2)*((1-sqrt(k1))/(1+sqrt(k1)));
  q= q0 + 2*q0^5 + 15*q0^9 + 150*q0^13; #(....)
  D=(10^(0.1*Rs)-1)/(10^(0.1*Rp)-1);

  n=ceil(log10(16*D)/log10(1/q));

endfunction

%% Test input validation
%!error ellipord ()
%!error ellipord (.1)
%!error ellipord (.1, .2)
%!error ellipord (.1, .2, 3)
%!error ellipord (.1, .2, 3, 4, 5)
%!error ellipord ([.1 .1], [.2 .2], 3, 4)
%!error ellipord ([.1 .2], [.5 .6], 3, 4)
%!error ellipord ([.1 .5], [.2 .6], 3, 4)
