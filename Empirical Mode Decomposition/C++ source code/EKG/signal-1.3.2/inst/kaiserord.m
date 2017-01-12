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
## @deftypefn  {Function File} {[@var{n}, @var{Wn}, @var{beta}, @var{ftype}] =} kaiserord (@var{f}, @var{m}, @var{dev})
## @deftypefnx {Function File} {[@dots{}] =} kaiserord (@var{f}, @var{m}, @var{dev}, @var{fs})
##
## Return the parameters needed to produce a filter of the desired
## specification from a Kaiser window.  The vector @var{f} contains pairs of
## frequency band edges in the range [0,1].  The vector @var{m} specifies the
## magnitude response for each band.  The values of @var{m} must be zero for
## all stop bands and must have the same magnitude for all pass bands. The
## deviation of the filter @var{dev} can be specified as a scalar or a vector
## of the same length as @var{m}.  The optional sampling rate @var{fs} can be
## used to indicate that @var{f} is in Hz in the range [0,@var{fs}/2].
##
## The returned value @var{n} is the required order of the filter (the length
## of the filter minus 1).  The vector @var{Wn} contains the band edges of
## the filter suitable for passing to @code{fir1}.  The value @var{beta} is
## the parameter of the Kaiser window of length @var{n}+1 to shape the filter.
## The string @var{ftype} contains the type of filter to specify to
## @code{fir1}.
##
## The Kaiser window parameters n and beta are computed from the
## relation between ripple (A=-20*log10(dev)) and transition width
## (dw in radians) discovered empirically by Kaiser:
##
## @example
## @group
##           / 0.1102(A-8.7)                        A > 50
##    beta = | 0.5842(A-21)^0.4 + 0.07886(A-21)     21 <= A <= 50
##           \ 0.0                                  A < 21
##
##    n = (A-8)/(2.285 dw)
## @end group
## @end example
##
## Example:
## @example
## @group
## [n, w, beta, ftype] = kaiserord ([1000, 1200], [1, 0], [0.05, 0.05], 11025);
## b = fir1 (n, w, kaiser (n+1, beta), ftype, "noscale");
## freqz (b, 1, [], 11025);
## @end group
## @end example
## @seealso{fir1, kaiser}
## @end deftypefn

## FIXME: order is underestimated for the final test case: 2 stop bands.

function [n, w, beta, ftype] = kaiserord(f, m, dev, fs)

  if (nargin<2 || nargin>4)
    print_usage;
  endif

  ## default sampling rate parameter
  if nargin<4, fs=2; endif

  ## parameter checking
  if length(f)!=2*length(m)-2
    error("kaiserord must have one magnitude for each frequency band");
  endif
  if any(m(1:length(m)-2)!=m(3:length(m)))
    error("kaiserord pass and stop bands must be strictly alternating");
  endif
  if length(dev)!=length(m) && length(dev)!=1
    error("kaiserord must have one deviation for each frequency band");
  endif
  dev = min(dev);
  if dev <= 0, error("kaiserord must have dev>0"); endif

  ## use midpoints of the transition region for band edges
  w = (f(1:2:length(f))+f(2:2:length(f)))/fs;

  ## determine ftype
  if length(w) == 1
    if m(1)>m(2), ftype='low'; else ftype='high'; endif
  elseif length(w) == 2
    if m(1)>m(2), ftype='stop'; else ftype='pass'; endif
  else
    if m(1)>m(2), ftype='DC-1'; else ftype='DC-0'; endif
  endif

  ## compute beta from dev
  A = -20*log10(dev);
  if (A > 50)
    beta = 0.1102*(A-8.7);
  elseif (A >= 21)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
  else
    beta = 0.0;
  endif

  ## compute n from beta and dev
  dw = 2*pi*min(f(2:2:length(f))-f(1:2:length(f)))/fs;
  n = max(1,ceil((A-8)/(2.285*dw)));

  ## if last band is high, make sure the order of the filter is even.
  if ((m(1)>m(2)) == (rem(length(w),2)==0)) && rem(n,2)==1, n = n+1; endif

endfunction

%!demo
%! Fs = 11025;
%! for i=1:4
%!   if i==1,
%!     subplot(221); bands=[1200, 1500]; mag=[1, 0]; dev=[0.1, 0.1];
%!   elseif i==2
%!     subplot(222); bands=[1000, 1500]; mag=[0, 1]; dev=[0.1, 0.1];
%!   elseif i==3
%!     subplot(223); bands=[1000, 1200, 3000, 3500]; mag=[0, 1, 0]; dev=0.1;
%!   elseif i==4
%!     subplot(224); bands=100*[10, 13, 15, 20, 30, 33, 35, 40];
%!     mag=[1, 0, 1, 0, 1]; dev=0.05;
%!   endif
%!   [n, w, beta, ftype] = kaiserord(bands, mag, dev, Fs);
%!   d=max(1,fix(n/10));
%!   if mag(length(mag))==1 && rem(d,2)==1, d=d+1; endif
%!   [h, f] = freqz(fir1(n,w,ftype,kaiser(n+1,beta),'noscale'),1,[],Fs);
%!   hm = freqz(fir1(n-d,w,ftype,kaiser(n-d+1,beta),'noscale'),1,[],Fs);
%!   plot(f,abs(hm),sprintf("r;order %d;",n-d), ...
%!        f,abs(h), sprintf("b;order %d;",n));
%!   b = [0, bands, Fs/2]; hold on;
%!   for i=2:2:length(b),
%!     hi=mag(i/2)+dev(1); lo=max(mag(i/2)-dev(1),0);
%!     plot([b(i-1), b(i), b(i), b(i-1), b(i-1)],[hi, hi, lo, lo, hi],"c;;");
%!   endfor; hold off;
%! endfor
%!
%! %--------------------------------------------------------------
%! % A filter meets the specifications if its frequency response
%! % passes through the ends of the criteria boxes, and fails if
%! % it passes through the top or the bottom.  The criteria are
%! % met precisely if the frequency response only passes through
%! % the corners of the boxes.  The blue line is the filter order
%! % returned by kaiserord, and the red line is some lower filter
%! % order.  Confirm that the blue filter meets the criteria and
%! % the red line fails.

## FIXME: extend demo to show detail at criteria box corners
