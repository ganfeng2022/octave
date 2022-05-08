## Copyright (C) Feng Gan <cesgf@mail.sysu.edu.cn; sysucesgf@163.com>
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
##

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{evolvingEigenValues}] =} fsmwefa (@var{X}, @var{winSize})
## Fixed size moving window evolving factor analysis
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X}                   --- a spectral matrix whose each row is a spectrum.
## @item
## @code{winSize}             --- size of the moving window.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{evolvingEigenValues} --- the eigen values in evolving process.
## @end itemize
## @end deftypefn
##
## References:
## - 1. Keller, H. R. and Massart, D. L. Anal. Chim. Acta, 1991, 246, 379 - 390.
##
## Lastest Revision: 2015-09-07
## Create date: 2008-07-08
function [evolvingEigenValues] = fsmwefa(X,winSize)
  evolvingEigenValues = [];
  if (nargin < 2)
    winSize = 5;  
  end
  [mRows, nCols] = size(X);
  if winSize > (min([mRows nCols]))/2
    error('You set a too big window size!');
  end
  for i = 1:(mRows - winSize + 1)
    subX = X(i:(i+winSize-1),:);
    subX = subX * subX';  
    eigenValues = svd(subX);
    eigenValues = log10(eigenValues); 
    evolvingEigenValues = [evolvingEigenValues; eigenValues'];
  end
endfunction
%!demo
%! load ./Data/masartdata.mat;
%! [ev] = fsmwefa(X,6);
%! plot(ev)





 
