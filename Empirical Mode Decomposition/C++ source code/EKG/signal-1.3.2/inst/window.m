## Copyright (C) 2008  David Bateman
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
## @deftypefn  {Function File} {@var{w} =} window (@var{f}, @var{m})
## @deftypefnx {Function File} {@var{w} =} window (@var{f}, @var{m}, @var{opts})
## Create an @var{m}-point window from the function @var{f}.  The function
## @var{f} can be for example @code{@@blackman}.  Any additional
## arguments @var{opt} are passed to the windowing function.
## @end deftypefn

function wout = window (f, m, varargin)

  if (nargin == 0)
    error ("window: UI tool not supported");
  elseif (nargin < 2)
    print_usage ();
  else
    w = feval (f, m, varargin{:});
    if (nargout > 0)
      wout = w;
    endif
  endif

endfunction

%!assert (window (@bartlett, 16), window ("bartlett", 16))
%!assert (window (@hamming, 16), window ("hamming", 16))
%!assert (window (@hanning, 16), window ("hanning", 16))
%!assert (window (@triang, 16), window ("triang", 16))

%% Test input validation
%!error window ()
%!error window (1)
%!error window ("hanning")
