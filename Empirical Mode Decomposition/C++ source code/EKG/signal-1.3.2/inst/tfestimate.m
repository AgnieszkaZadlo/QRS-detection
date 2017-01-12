## Copyright (C) 2006 Peter V. Lanspeary <pvl@mecheng.adelaide.edu.au>
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
## @deftypefn  {Function File} {} tfestimate (@var{x}, @var{y})
## @deftypefnx {Function File} {} tfestimate (@var{x}, @var{y}, @var{window})
## @deftypefnx {Function File} {} tfestimate (@var{x}, @var{y}, @var{window}, @var{overlap})
## @deftypefnx {Function File} {} tfestimate (@var{x}, @var{y}, @var{window}, @var{overlap}, @var{Nfft})
## @deftypefnx {Function File} {} tfestimate (@var{x}, @var{y}, @var{window}, @var{overlap}, @var{Nfft}, @var{Fs})
## @deftypefnx {Function File} {} tfestimate (@var{x}, @var{y}, @var{window}, @var{overlap}, @var{Nfft}, @var{Fs}, @var{range})
## @deftypefnx {Function File} {[@var{Pxx}, @var{freq}] =} tfestimate (@dots{})
##
## Estimate transfer function of system with input @var{x} and output @var{y}.
## Use the Welch (1967) periodogram/FFT method.
## @seealso{pwelch}
## @end deftypefn

function varargout = tfestimate(varargin)

  ##
  ## Check fixed argument
  if (nargin < 2 || nargin > 7)
    print_usage ();
  endif
  nvarargin = length(varargin);
  ## remove any pwelch RESULT args and add 'cross'
  for iarg=1:nvarargin
    arg = varargin{iarg};
    if ( ~isempty(arg) && ischar(arg) && ( strcmp(arg,'power') || ...
           strcmp(arg,'cross') || strcmp(arg,'trans') || ...
           strcmp(arg,'coher') || strcmp(arg,'ypower') ))
      varargin{iarg} = [];
    endif
  endfor
  varargin{nvarargin+1} = 'trans';
  ##
  if ( nargout==0 )
    pwelch(varargin{:});
  elseif ( nargout==1 )
    Pxx = pwelch(varargin{:});
    varargout{1} = Pxx;
  elseif ( nargout>=2 )
    [Pxx,f] = pwelch(varargin{:});
    varargout{1} = Pxx;
    varargout{2} = f;
  endif

endfunction
