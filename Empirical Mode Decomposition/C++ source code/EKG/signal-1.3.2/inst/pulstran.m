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
## @deftypefn  {Function File} {@var{y} =} pulstran (@var{t}, @var{d}, @var{func}, @dots{})
## @deftypefnx {Function File} {@var{y} =} pulstran (@var{t}, @var{d}, @var{p})
## @deftypefnx {Function File} {@var{y} =} pulstran (@var{t}, @var{d}, @var{p}, @var{Fs})
## @deftypefnx {Function File} {@var{y} =} pulstran (@var{t}, @var{d}, @var{p}, @var{Fs}, @var{method})
##
## Generate the signal y=sum(func(t+d,...)) for each d.  If d is a
## matrix of two columns, the first column is the delay d and the second
## column is the amplitude a, and y=sum(a*func(t+d)) for each d,a.
## Clearly, func must be a function which accepts a vector of times.
## Any extra arguments needed for the function must be tagged on the end.
##
## Example:
## @example
## fs = 11025;  # arbitrary sample rate
## f0 = 100;    # pulse train sample rate
## w = 0.001;   # pulse width of 1 millisecond
## auplot (pulstran (0:1/fs:0.1, 0:1/f0:0.1, "rectpuls", w), fs);
## @end example
##
## If instead of a function name you supply a pulse shape sampled at
## frequency Fs (default 1 Hz),  an interpolated version of the pulse
## is added at each delay d.  The interpolation stays within the the
## time range of the delayed pulse.  The interpolation method defaults
## to linear, but it can be any interpolation method accepted by the
## function interp1.
##
## Example:
## @example
## fs = 11025;      # arbitrary sample rate
## f0 = 100;        # pulse train sample rate
## w = boxcar(10);  # pulse width of 1 millisecond at 10 kHz
## auplot (pulstran (0:1/fs:0.1, 0:1/f0:0.1, w, 10000), fs);
## @end example
## @end deftypefn

## FIXME: Make it faster.  It is currently unusable for anything real.
##        It may not be possible to speed it up with the present interface.
##        See speech/voice.m for a better way.

## Note that pulstran can be used for some pretty strange things such
## as simple band-limited interpolation:
##     xf = 0:0.05:10; yf = sin(2*pi*xf/5);
##     xp = 0:10; yp = sin(2*pi*xp/5); # .2 Hz sine sampled every second
##     s = pulstran(xf, [xp, yp],'sinc');
##     plot(f, yf, ";original;", xf, s, ";sinc;",xp,yp,"*;;");
## You wouldn't want to do this in practice since it is expensive, and
## since it works much better with a windowed sinc function, at least
## for short samples.

function y = pulstran(t, d, pulse, varargin)

  if nargin<3 || (!ischar(pulse) && nargin>5)
    print_usage;
  endif
  y = zeros(size(t));
  if isempty(y), return; endif
  if rows(d) == 1, d=d'; endif
  if columns(d) == 2,
    a=d(:,2);
  else
    a=ones(rows(d),1);
  endif
  if ischar(pulse)
    ## apply function t+d for all d
    for i=1:rows(d)
      y = y+a(i)*feval(pulse,t-d(i,1),varargin{:});
    endfor
  else
    ## interpolate each pulse at the specified times
    Fs = 1; method = 'linear';
    if nargin==4
      arg=varargin{1};
      if ischar(arg),
        method=arg;
      else
        Fs = arg;
      endif
    elseif nargin==5
      Fs = varargin{1};
      method = varargin{2};
    endif
    span = (length(pulse)-1)/Fs;
    t_pulse = (0:length(pulse)-1)/Fs;
    for i=1:rows(d)
      dt = t-d(i,1);
      idx = find(dt>=0 & dt<=span);
      y(idx) = y(idx) + a(i)*interp1(t_pulse, pulse, dt(idx), method);
    endfor
  endif

endfunction

%!error pulstran
%!error pulstran(1,2,3,4,5,6)

%!## parameter size and shape checking
%!shared t,d
%! t = 0:0.01:1; d=0:0.1:1;
%!assert (isempty(pulstran([], d, 'sin')));
%!assert (pulstran(t, [], 'sin'), zeros(size(t)));
%!assert (isempty(pulstran([], d, boxcar(5))));
%!assert (pulstran(t, [], boxcar(5)), zeros(size(t)));
%!assert (size(pulstran(t,d,'sin')), size(t));
%!assert (size(pulstran(t,d','sin')), size(t));
%!assert (size(pulstran(t',d,'sin')), size(t'));
%!assert (size(pulstran(t,d','sin')), size(t));

%!demo
%! fs = 11025;                   # arbitrary sample rate
%! f0 = 100;                     # pulse train sample rate
%! w = 0.003;                    # pulse width of 3 milliseconds
%! t = 0:1/fs:0.1; d=0:1/f0:0.1; # define sample times and pulse times
%! a = hanning(length(d));       # define pulse amplitudes
%!
%! subplot(221);
%! x = pulstran(t', d', 'rectpuls', w);
%! plot([0:length(x)-1]*1000/fs, x);
%! hold on; plot(d*1000,ones(size(d)),'g*;pulse;'); hold off;
%! ylabel("amplitude"); xlabel("time (ms)");
%! title("rectpuls");
%!
%! subplot(223);
%! x = pulstran(f0*t, [f0*d', a], 'sinc');
%! plot([0:length(x)-1]*1000/fs, x);
%! hold on; plot(d*1000,a,'g*;pulse;'); hold off;
%! ylabel("amplitude"); xlabel("time (ms)");
%! title("sinc => band limited interpolation");
%!
%! subplot(222);
%! pulse = boxcar(30);  # pulse width of 3 ms at 10 kHz
%! x = pulstran(t, d', pulse, 10000);
%! plot([0:length(x)-1]*1000/fs, x);
%! hold on; plot(d*1000,ones(size(d)),'g*;pulse;'); hold off;
%! ylabel("amplitude"); xlabel("time (ms)");
%! title("interpolated boxcar");
%!
%! subplot(224);
%! pulse = sin(2*pi*[0:0.0001:w]/w).*[w:-0.0001:0];
%! x = pulstran(t', [d', a], pulse', 10000);
%! plot([0:length(x)-1]*1000/fs, x);
%! hold on; plot(d*1000,a*w,'g*;pulse;'); hold off; title("");
%! ylabel("amplitude"); xlabel("time (ms)");
%! title("interpolated asymmetric sin");
%!
%! %----------------------------------------------------------
%! % Should see (1) rectangular pulses centered on *,
%! %            (2) rectangular pulses to the right of *,
%! %            (3) smooth interpolation between the *'s, and
%! %            (4) asymmetric sines to the right of *
