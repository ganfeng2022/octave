
## Copyright (C) <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
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
## @deftypefn {Function File} {[@var{residuals}]} = plsNipalsCv ([@var{c},@var{Y},@var{nPrinComp},@var{itMax},@var{tol}])
## Cross validation
##
## Input arguments:
##
## @itemize
## @code{c}           --- property vector.
## @item
## @code{Y}           --- spectral matrix whose each row is a sample.
## @item
## @code{nPrinComp}   --- number principal component.
## @item
## @code{itMax}       --- maximum iteration.
## @item
## @code{tol}         --- convergence criterion.
## @end itemize
## 
## Return values
##
## @itemize
## @item
## @code{residuals}   --- residuals at different nPrinComp.
## @end itemize
##
## @seealso{pls1_nipals, pls1_pred.}
## @end deftypefn

## Author:  Gan, F.
## create date:     2017-10-21
## latest revision: 2017-10-21

function [residuals] = plsNipalsCv(c,Y,nPrinComp,itMax,tol)

  if (nargin < 5)
    error('Please see demo.');
  endif
  if (size(c,2) > 1)
    error('This is for PLS1.');
  endif
  [mRows,nCols] = size(Y);
  nPrinComp = min([mRows,nCols,nPrinComp]);
  residuals = zeros(nPrinComp,1);
  for i = 1:mRows
    ## cross one line
    range = [1:(i-1),(i+1):mRows];
    subY = Y(range,:); 
    subc = c(range,1);
    [meancentY,barY] = meancent(subY,'col');
    [meancentc,barc] = meancent(subc,'col');
    [W,P,T,U,b,q] = plsNipals(meancentc,meancentY,nPrinComp,itMax,tol);
    for j = 1:nPrinComp
      [hatc] = plsNipalsPred(Y(i,:),W,P,b,q,j,barY,barc);
      res = hatc - c(i);
      residuals(j) = residuals(j) + res * res; 
    endfor
  endfor
  residuals = residuals ./ mRows;
endfunction

%!demo
%! load ./Data/CornModelSp.dat;
%! load ./Data/CornModelProp.dat;
%! pcn = 10;
%! itmax = 100;
%! tol = 1e-6;
%! [cv] = plsNipalsCv(CornModelProp,CornModelSp,pcn,itmax,tol);
%! figure(1),clf('reset');
%! plot(cv);hold on;plot(cv,'o');
%! xlabel('Principal component number');
%! ylabel('Residuals');

