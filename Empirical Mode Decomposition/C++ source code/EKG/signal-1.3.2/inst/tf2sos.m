## Copyright (C) 2005 Julius O. Smith III <jos@ccrma.stanford.edu>
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
## @deftypefn  {Function File} {[@var{sos}, @var{g}] =} tf2sos (@var{b}, @var{a})
## @deftypefnx {Function File} {@var{sos} =} tf2sos (@var{b}, @var{a})
## Convert direct-form filter coefficients to series second-order sections.
##
## INPUTS:
##
## @var{b} and @var{a} are vectors specifying the digital filter
## @math{H(z) = B(z)/A(z)}.  See @code{filter} for documentation of the @var{b}
## and @var{a} filter coefficients.
##
## RETURNED:
## @itemize
## @item
## @var{sos} = matrix of series second-order sections, one per row:
## @example
## @var{sos} = [@var{b1}.' @var{a1}.'; ...; @var{bn}.' @var{an}.']
## @end example
## where
## @code{@var{B1}.' = [b0 b1 b2] and @var{A1}.' = [1 a1 a2]} for
## section 1, etc.  The b0 entry must be nonzero for each section (zeros at
## infinity not supported).
## @item
## @var{g} is an overall gain factor that effectively scales
## any one of the @var{B}i vectors.
## @end itemize
##
## If called with only one output argument, the overall filter gain is
## applied to the first second-order section in the matrix @var{sos}.
##
## EXAMPLE:
## @example
## B = [1 0 0 0 0 1];
## A = [1 0 0 0 0 .9];
## [sos, g] = tf2sos (B, A)
##
## sos =
##
##    1.00000   0.61803   1.00000   1.00000   0.60515   0.95873
##    1.00000  -1.61803   1.00000   1.00000  -1.58430   0.95873
##    1.00000   1.00000  -0.00000   1.00000   0.97915  -0.00000
##
## g = 1
## @end example
##
## @seealso{sos2tf, zp2sos, sos2pz, zp2tf, tf2zp}
## @end deftypefn

function [sos,g] = tf2sos (B, A)

  [z,p,k] = tf2zp(B(:)',A(:)');
  if (nargout < 2)
    sos = zp2sos(z,p,k);
  else
    [sos,g] = zp2sos(z,p,k);
  endif

endfunction

%!test
%! B=[1 0 0 0 0 1]; A=[1 0 0 0 0 .9];
%! [sos,g] = tf2sos(B,A);
%! [Bh,Ah] = sos2tf(sos,g);
%! assert({Bh,Ah},{B,A},100*eps);
