## Copyright (c) 2013 Rob Sykes <robs@users.sourceforge.net>
##
## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
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
## @deftypefn  {Function File} {[@var{w}, @var{xmu}] =} ultrwin (@var{m}, @var{mu}, @var{beta})
## @deftypefnx {Function File} {[@var{w}, @var{xmu}] =} ultrwin (@var{m}, @var{mu}, @var{att}, "att")
## @deftypefnx {Function File} {[@var{w}, @var{xmu}] =} ultrwin (@var{m}, @var{mu}, @var{latt}, "latt")
## @deftypefnx {Function File} {@var{w} =} ultrwin (@var{m}, @var{mu}, @var{xmu}, "xmu")
## Return the coefficients of an Ultraspherical window of length @var{m}.
## The parameter @var{mu} controls the window's Fourier transform's side-lobe
## to side-lobe ratio, and the third given parameter controls the transform's
## main-lobe width/side-lobe-ratio; normalize @var{w} such that the central
## coefficient(s) value is unitary.
##
## By default, the third parameter is @var{beta}, which sets the main lobe width
## to @var{beta} times that of a rectangular window.  Alternatively, giving
## @var{att} or @var{latt} sets the ripple ratio at the first or last side-lobe
## respectively, or giving @var{xmu} sets the (un-normalized) window's Fourier
## transform according to its canonical definition:
##
## @verbatim
##              (MU)
##      W(k) = C   [ XMU cos(pi k/M) ],  k = 0, 1, ..., M-1,
##              M-1
## @end verbatim
##
## where C is the Ultraspherical (a.k.a. Gegenbauer) polynomial, which can be
## defined using the recurrence relationship:
##
## @verbatim
##       (l)    1                  (l)                    (l)
##      C (x) = - [ 2x(m + l - 1) C   (x) - (m + 2l - 2) C   (x) ]
##       m      m                  m-1                    m-2
##
##                                 (l)        (l)
##      for m an integer > 1, and C (x) = 1, C (x) = 2lx.
##                                 0          1
## @end verbatim
##
## For given @var{beta}, @var{att}, or @var{latt}, the corresponding
## (determined) value of @var{xmu} is also returned.
##
## The Dolph-Chebyshev and Saramaki windows are special cases of the
## Ultraspherical window, with @var{mu} set to 0 and 1 respectively.  Note that
## when not giving @var{xmu}, stability issues may occur with @var{mu} <= -1.5.
## For further information about the window, see
##
## @itemize @bullet
## @item
## Kabal, P., 2009: Time Windows for Linear Prediction of Speech.
## Technical Report, Dept. Elec. & Comp. Eng., McGill University.
## @item
## Bergen, S., Antoniou, A., 2004: Design of Ultraspherical Window
## Functions with Prescribed Spectral Characteristics. Proc. JASP, 13/13,
## pp. 2053-2065.
## @item
## Streit, R., 1984: A two-parameter family of weights for nonrecursive
## digital filters and antennas. Trans. ASSP, 32, pp. 108-118.
## @end itemize
## @seealso{chebwin, kaiser}
## @end deftypefn

function [w, xmu] = ultrwin (m, mu, par, key = "beta", norm = 0)
  ## This list of parameter types must be kept in sync with the enum order.
  types = {"xmu", "beta", "att", "latt"};
  type = [];
  if (ischar (key))
    type = find (strncmpi (key, types, numel (key)));
  endif

  if (nargin < 3 || nargin > 5)
    print_usage ();
  elseif (! (isscalar (m) && (m == fix (m)) && (m > 0)))
    error ("ultrwin: M must be a positive integer");
  elseif (! (isscalar (mu) && isreal (mu)))
    error ("ultrwin: MU must be a real scalar");
  elseif (! ischar (key))
    error ("ultrwin: parameter type must be a string");
  elseif (isempty (type))
    error ("ultrwin: invalid parameter type '%s'", key);
  elseif (! (isscalar (par) && isreal (par)))
    error (["ultrwin: ", upper (types(type)), " must be a real scalar"]);
  elseif (! (isscalar (norm) && norm == fix (norm) && norm >= 0)) # Alt. norms; WIP
    error ("ultrwin: NORM must be a non-negative integer");
  endif

  [w, xmu] = __ultrwin__(m, mu, par, type-1, norm);

endfunction

%!test
%! assert(ultrwin(100, 1, 1), ones(100, 1), 1e-14);

%!test
%! L = 201; xmu = 1.01; m = L-1;
%! for mu = -1.35:.3:1.35
%!   x = xmu*cos([0:m]*pi/L);
%!   C(2,:) = 2*mu*x; C(1,:) = 1;
%!   for k = 2:m; C(k+1,:) = 2*(k+mu-1)/k*x.*C(k,:) - (k+2*mu-2)/k*C(k-1,:); end
%!   b = real(ifft(C(m+1,:))); b = b(m/2+2:L)/b(1);
%!   assert(ultrwin(L, mu, xmu, "x")', [b 1 fliplr(b)], 1e-12);
%! end

%!test
%! b = [
%!   5.7962919401511820e-03
%!   1.6086991349967078e-02
%!   3.6019014684117417e-02
%!   6.8897525451558125e-02
%!   1.1802364384553447e-01
%!   1.8566749737411145e-01
%!   2.7234740630826737e-01
%!   3.7625460141456091e-01
%!   4.9297108901880221e-01
%!   6.1558961695849457e-01
%!   7.3527571856983598e-01
%!   8.4222550739092694e-01
%!   9.2688779484512085e-01
%!   9.8125497127708561e-01]';
%! [w xmu] = ultrwin(29, 0, 3);
%! assert(w', [b 1 fliplr(b)], 1e-14);
%! assert(xmu, 1.053578297819277, 1e-14);

%!test
%! b = [
%!   2.9953636903962466e-02
%!   7.6096450051659603e-02
%!   1.5207129867916891e-01
%!   2.5906995366355179e-01
%!   3.9341065451220536e-01
%!   5.4533014012036929e-01
%!   6.9975915071207051e-01
%!   8.3851052636906720e-01
%!   9.4345733548690369e-01]';
%! assert(ultrwin(20, .5, 50, "a")', [b 1 1 fliplr(b)], 1e-14);

%!test
%! b = [
%!   1.0159906492322712e-01
%!   1.4456358609406283e-01
%!   2.4781689516201011e-01
%!   3.7237015168857646e-01
%!   5.1296973026690407e-01
%!   6.5799041448113671e-01
%!   7.9299087042967320e-01
%!   9.0299778924260576e-01
%!   9.7496213649820296e-01]';
%! assert(ultrwin(19, -.4, 40, "l")', [b 1 fliplr(b)], 1e-14);

%!demo
%! w=ultrwin(120, -1, 40, "l"); [W,f]=freqz(w); clf
%! subplot(2,1,1); plot(f/pi, 20*log10(W/abs(W(1)))); grid; axis([0 1 -90 0])
%! subplot(2,1,2); plot(0:length(w)-1, w); grid
%! %-----------------------------------------------------------
%! % Figure shows an Ultraspherical window with MU=-1, LATT=40:
%! % frequency domain above, time domain below.

%!demo
%! c="krbm"; clf; subplot(2, 1, 1)
%! for beta=2:5
%!   w=ultrwin(80, -.5, beta); [W,f]=freqz(w);
%!   plot(f/pi, 20*log10(W/abs(W(1))), c(1+mod(beta, length(c)))); hold on
%! end; grid; axis([0 1 -140 0]); hold off
%! subplot(2, 1, 2);
%! for n=2:10
%!   w=ultrwin(n*20, 1, 3); [W,f]=freqz(w,1,2^11);
%!   plot(f/pi, 20*log10(W/abs(W(1))), c(1+mod(n, length(c)))); hold on
%! end; grid; axis([0 .2 -100 0]); hold off
%! %--------------------------------------------------
%! % Figure shows transfers of Ultraspherical windows:
%! % above: varying BETA with fixed N & MU,
%! % below: varying N with fixed MU & BETA.

%!demo
%! c="krbm"; clf; subplot(2, 1, 1)
%! for j=0:4
%!   w=ultrwin(80, j*.6-1.2, 50, "a"); [W,f]=freqz(w);
%!   plot(f/pi, 20*log10(W/abs(W(1))), c(1+mod(j, length(c)))); hold on
%! end; grid; axis([0 1 -100 0]); hold off
%! subplot(2, 1, 2);
%! for j=4:-1:0
%!   w=ultrwin(80, j*.75-1.5, 50, "l"); [W,f]=freqz(w);
%!   plot(f/pi, 20*log10(W/abs(W(1))), c(1+mod(j, length(c)))); hold on
%! end; grid; axis([0 1 -100 0]); hold off
%! %--------------------------------------------------
%! % Figure shows transfers of Ultraspherical windows:
%! % above: varying MU with fixed N & ATT,
%! % below: varying MU with fixed N & LATT.

%!demo
%! clf; a=[.8 2 -115 5]; fc=1.1/pi; l="labelxy";
%! for k=1:3; switch (k); case 1; w=kaiser(L=159, 7.91);
%!   case 2; w=ultrwin(L=165, 0, 2.73); case 3; w=ultrwin(L=153, .5, 2.6); end
%!   subplot(3, 1, 4-k); f=[1:(L-1)/2]*pi;f=sin(fc*f)./f; f=[fliplr(f) fc f]';
%!   [h,f]=freqz(w.*f,1,2^14); plot(f,20*log10(h)); grid; axis(a,l); l="labely";
%! end
%! %-----------------------------------------------------------
%! % Figure shows example lowpass filter design (Fp=1, Fs=1.2
%! % rad/s, att=80 dB) and comparison with other windows.  From
%! % top to bottom: Ultraspherical, Dolph-Chebyshev, and Kaiser
%! % windows, with lengths 153, 165, and 159 respectively.
