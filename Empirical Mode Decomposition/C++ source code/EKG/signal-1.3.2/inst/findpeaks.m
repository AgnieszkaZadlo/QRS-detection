## Copyright (c) 2012 Juan Pablo Carbajal <carbajal@ifi.uzh.ch>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{pks}, @var{loc}, @var{extra}] =} findpeaks (@var{data})
## @deftypefnx {Function File} {@dots{} =} findpeaks (@dots{}, @var{property}, @var{value})
## @deftypefnx {Function File} {@dots{} =} findpeaks (@dots{}, @asis{"DoubleSided"})
## Finds peaks on @var{data}.
##
## Peaks of a positive array of data are defined as local maxima. For
## double-sided data, they are maxima of the positive part and minima of
## the negative part. @var{data} is expected to be a single column
## vector.
##
## The function returns the value of @var{data} at the peaks in
## @var{pks}. The index indicating their position is returned in
## @var{loc}.
##
## The third output argument is a structure with additional information:
##
## @table @asis
## @item "parabol"
## A structure containing the parabola fitted to each returned peak. The
## structure has two fields, @asis{"x"} and @asis{"pp"}. The field
## @asis{"pp"} contains the coefficients of the 2nd degree polynomial
## and @asis{"x"} the extrema of the intercal here it was fitted.
##
## @item "height"
## The estimated height of the returned peaks (in units of @var{data}).
##
## @item "baseline"
## The height at which the roots of the returned peaks were calculated
## (in units of @var{data}).
##
## @item "roots"
## The abscissa values (in index units) at which the parabola fitted to
## each of the returned peaks crosses the @asis{"baseline"} value. The
## width of the peak is calculated by @command{diff(roots)}.
## @end table
##
## This function accepts property-value pair given in the list below:
##
## @table @asis
##
## @item "MinPeakHeight"
## Minimum peak height (positive scalar). Only peaks that exceed this
## value will be returned. For data taking positive and negative values
## use the option "DoubleSided". Default value @code{2*std (abs (detrend
## (data,0)))}.
##
## @item "MinPeakDistance"
## Minimum separation between (positive integer). Peaks separated by
## less than this distance are considered a single peak. This distance
## is also used to fit a second order polynomial to the peaks to
## estimate their width, therefore it acts as a smoothing parameter.
## Default value 4.
##
## @item "MinPeakWidth"
## Minimum width of peaks (positive integer). The width of the peaks is
## estimated using a parabola fitted to the neighborhood of each peak.
## The neighborhood size is equal to the value of
## @asis{"MinPeakDistance"}. The width is evaluated at the half height
## of the peak with baseline at "MinPeakHeight". Default value 2.
##
## @item "DoubleSided"
## Tells the function that data takes positive and negative values. The
## base-line for the peaks is taken as the mean value of the function.
## This is equivalent as passing the absolute value of the data after
## removing the mean.
## @end table
##
## Run @command{demo findpeaks} to see some examples.
## @end deftypefn

function [pks idx varargout] = findpeaks (data, varargin)

  if (nargin < 1)
    print_usage ();
  endif

  if (! (isvector (data) && numel (data) >= 3))
    error ("findpeaks:InvalidArgument",
           "findpeaks: DATA must be a vector of at least 3 elements");
  endif

  transpose = (rows (data) == 1);

  if (transpose)
    data = data.';
  endif

  ## --- Parse arguments --- #
  __data__ = abs (detrend (data, 0));

  posscal = @(x) isscalar (x) && x >= 0;

  parser = inputParser ();
  parser.FunctionName = "findpeaks";

  ## FIXME: inputParser was first implemented in the general package in the
  ##        old @class type.  This allowed for a very similar interface to
  ##        Matlab but not quite equal.  classdef was then implemented in
  ##        Octave 4.0 release, which enabled inputParser to be implemented
  ##        properly.  However, this causes problem because we don't know
  ##        what implementation may be running.  A new version of the general
  ##        package is being released to avoid the two implementations to
  ##        co-exist.
  ##
  ##        To keep supporting older octave versions, we have an alternative
  ##        path that avoids inputParser.  And if inputParser is available,
  ##        we check what implementation is is, and act accordingly.

  ## Note that in Octave 4.0, inputParser is classdef and Octave behaves
  ## weird for it. which ("inputParser") will return empty (thinks its a
  ## builtin function).
  if (exist ("inputParser") == 2
      && isempty (strfind (which ("inputParser"),
                           ["@inputParser" filesep "inputParser.m"])))
    ## making use of classdef's inputParser ..
    parser.addParamValue ("MinPeakHeight", 2*std (__data__),posscal);
    parser.addParamValue ("MinPeakDistance", 4, posscal);
    parser.addParamValue ("MinPeakWidth", 2, posscal);
    parser.addSwitch ("DoubleSided");
    parser.parse (varargin{:});
    minH      = parser.Results.MinPeakHeight;
    minD      = parser.Results.MinPeakDistance;
    minW      = parser.Results.MinPeakWidth;
    dSided    = parser.Results.DoubleSided;
  else
    ## either old @inputParser or no inputParser at all...
    lvarargin = lower (varargin);

    ds = strcmpi (lvarargin, "DoubleSided");
    if (any (ds))
      dSided = true;
      lvarargin(ds) = [];
    else
      dSided = false;
    endif

    [~, minH, minD, minW] = parseparams (lvarargin,
                                         "minpeakheight", 2 * std (__data__),
                                         "minpeakdistance", 4,
                                         "minpeakwidth", 2);
    if (! posscal (minH))
      error ("findpeaks: MinPeakHeight must be a positive scalar");
    elseif (! posscal (minD))
      error ("findpeaks: MinPeakDistance must be a positive scalar");
    elseif (! posscal (minW))
      error ("findpeaks: MinPeakWidth must be a positive scalar");
    endif
  endif


  if (dSided)
    [data, __data__] = deal (__data__, data);
  elseif (min (data) < 0)
    error ("findpeaks:InvalidArgument",
           'Data contains negative values. You may want to "DoubleSided" option');
  endif

  ## Rough estimates of first and second derivative
  df1 = diff (data, 1)([1; (1:end)']);
  df2 = diff (data, 2)([1; 1; (1:end)']);

  ## check for changes of sign of 1st derivative and negativity of 2nd
  ## derivative.
  ## <= in 1st derivative includes the case of oversampled signals.
  idx = find (df1.*[df1(2:end); 0] <= 0 & [df2(2:end); 0] < 0);

  ## Get peaks that are beyond given height
  tf  = data(idx) > minH;
  idx = idx(tf);

  ## sort according to magnitude
  [~, tmp] = sort (data(idx), "descend");
  idx_s = idx(tmp);

  ## Treat peaks separated less than minD as one
  D = abs (bsxfun (@minus, idx_s, idx_s'));
  if (any (D(:) < minD))

    i = 1;
    peak = cell ();
    node2visit = 1:size(D,1);
    visited = [];
    idx_pruned = idx_s;

    ## debug
##    h = plot(1:length(data),data,"-",idx_s,data(idx_s),'.r',idx_s,data(idx_s),'.g');
##    set(h(3),"visible","off");

    while (! isempty (node2visit))

      d = D(node2visit(1),:);

      visited = [visited node2visit(1)];
      node2visit(1) = [];

      neighs  = setdiff (find (d < minD), visited);
      if (! isempty (neighs))
        ## debug
##        set(h(3),"xdata",idx_s(neighs),"ydata",data(idx_s(neighs)),"visible","on")
##        pause(0.2)
##        set(h(3),"visible","off");

        idx_pruned = setdiff (idx_pruned, idx_s(neighs));

        visited    = [visited neighs];
        node2visit = setdiff (node2visit, visited);

        ## debug
##        set(h(2),"xdata",idx_pruned,"ydata",data(idx_pruned))
##        pause
      endif

    endwhile
    idx = idx_pruned;
  endif

  extra = struct ("parabol", [], "height", [], "baseline", [], "roots", []);

  ## Estimate widths of peaks and filter for:
  ## width smaller than given.
  ## wrong concavity.
  ## not high enough
  ## data at peak is lower than parabola by 1%
  if (minW > 0)
    ## debug
##    h = plot(1:length(data),data,"-",idx,data(idx),'.r',...
##          idx,data(idx),'og',idx,data(idx),'-m');
##    set(h(4),"linewidth",2)
##    set(h(3:4),"visible","off");

    idx_pruned = idx;
    n  = numel (idx);
    np = numel (data);
    struct_count = 0;
    for i=1:n
      ind = (round (max(idx(i)-minD/2,1)) : ...
             round (min(idx(i)+minD/2,np)))';
      pp = polyfit (ind, data(ind), 2);
      H  = pp(3) - pp(2)^2/(4*pp(1));

      ## debug
##      x = linspace(ind(1)-1,ind(end)+1,10);
##      set(h(4),"xdata",x,"ydata",polyval(pp,x),"visible","on")
##      set(h(3),"xdata",ind,"ydata",data(ind),"visible","on")
##      pause(0.2)
##      set(h(3:4),"visible","off");

      rz = roots ([pp(1:2) pp(3)-mean([H,minH])]);
      width = abs (diff (rz));
      if (width < minW || pp(1) > 0 || H < minH || data(idx(i)) < 0.99*H)
        idx_pruned = setdiff (idx_pruned, idx(i));
      elseif (nargout > 2)
        struct_count++;
        extra.parabol(struct_count).x  = ind([1 end]);
        extra.parabol(struct_count).pp = pp;

        extra.roots(struct_count,1:2)= rz;
        extra.height(struct_count)   = H;
        extra.baseline(struct_count) = mean ([H minH]);
      endif

      ## debug
##      set(h(2),"xdata",idx_pruned,"ydata",data(idx_pruned))
##      pause(0.2)

    endfor
    idx = idx_pruned;
  endif


  if (dSided)
    pks = __data__(idx);
  else
    pks = data(idx);
  endif

  if (transpose)
    pks = pks.';
    idx = idx.';
  endif

  if (nargout() > 2)
    varargout{1} = extra;
  endif

endfunction

%!demo
%! t = 2*pi*linspace(0,1,1024)';
%! y = sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3);
%!
%! data1 = abs(y); # Positive values
%! [pks idx] = findpeaks(data1);
%!
%! data2 = y; # Double-sided
%! [pks2 idx2] = findpeaks(data2,"DoubleSided");
%! [pks3 idx3] = findpeaks(data2,"DoubleSided","MinPeakHeight",0.5);
%!
%! subplot(1,2,1)
%! plot(t,data1,t(idx),data1(idx),'.m')
%! subplot(1,2,2)
%! plot(t,data2,t(idx2),data2(idx2),".m;>2*std;",t(idx3),data2(idx3),"or;>0.1;")
%! legend("Location","NorthOutside","Orientation","horizontal")
%!
%! #----------------------------------------------------------------------------
%! # Finding the peaks of smooth data is not a big deal!

%!demo
%! t = 2*pi*linspace(0,1,1024)';
%! y = sin(3.14*t) + 0.5*cos(6.09*t) + 0.1*sin(10.11*t+1/6) + 0.1*sin(15.3*t+1/3);
%!
%! data = abs(y + 0.1*randn(length(y),1)); # Positive values + noise
%! [pks idx] = findpeaks(data,"MinPeakHeight",1);
%!
%! dt = t(2)-t(1);
%! [pks2 idx2] = findpeaks(data,"MinPeakHeight",1,...
%!                              "MinPeakDistance",round(0.5/dt));
%!
%! subplot(1,2,1)
%! plot(t,data,t(idx),data(idx),'.r')
%! subplot(1,2,2)
%! plot(t,data,t(idx2),data(idx2),'.r')
%!
%! #----------------------------------------------------------------------------
%! # Noisy data may need tuning of the parameters. In the 2nd example,
%! # MinPeakDistance is used as a smoother of the peaks.

%!assert (isempty (findpeaks ([1, 1, 1])))
%!assert (isempty (findpeaks ([1; 1; 1])))

## Test for bug #45056
%!test
%! ## Test input vector is an oversampled sinusoid with clipped peaks
%! x = min (3, cos (2*pi*[0:8000] ./ 600) + 2.01);
%! assert (! isempty (findpeaks (x)))

%% Test input validation
%!error findpeaks ()
%!error findpeaks (1)
%!error findpeaks ([1, 2])

## Failing test because we are not Matlab compatible
%!xtest assert (findpeaks ([34 134 353 64 134 14 56 67 234 143 64 575 8657]),
%!              [353 134 234])

