
## Copyright (C) Feng Gan <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{z}]} = whittakersmoother ([@var{x},@var{w},@var{lambda},@var{order}])
## A perfect smooter based on Whittaker's theory.
##
## Reference:
## Eilers, P. H. C. A perfect smoother. Anal. Chem. 2003. 75. 3631-3636. 
##
## Input arguments:
##
## @itemize
## @item
## @code{unsmoothSignal}  --- unsmooth singnal vector.
## @item
## @code{singnalWeight}   --- signal weight vector.
## @item
## @code{lambda}          --- smoothing parameter, the bigger the smoother.
## @item
## @code{order}           --- order of difference.
## @end itemize
## 
## Return values
##
## @itemize
## @item
## @code{smoothedSignal}  --- smoothed signal vector.
## @end itemize
##
## @seealso{functions}
## @end deftypefn

## Author:  Feng GAN
## create date:     2017-09-01
## latest revision: 2017-09-01

function [smoothedSignal] = whittakersmoother(unsmoothSignal,signalWeight,lambda,order)

  if nargin < 4
    error('Please see demo part.');
  endif

  mRows = length(unsmoothSignal);
  E = speye(mRows);
  D = diff(E, order);
  W = spdiags(signalWeight,0,mRows,mRows);
  C = chol(W + lambda * D' * D);
  smoothedSignal = C \ (C' \ (signalWeight .* unsmoothSignal));
  
endfunction

%!demo
%! load ./Data/noisedata.txt;
%! x = noisedata(:,2);
%! w = ones(size(x));
%! lambda = 1e6;
%! order = 2;
%! [z] = whittakersmoother(x,w,lambda,order);
%! figure(1),clf('reset'),plot(z)

