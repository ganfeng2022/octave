
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
## @deftypefn {Function File} {[@var{W},@var{P},@var{T},@var{U},@var{b},@var{q}]} = plsNipals ([@var{y},@var{X},@var{pcn},@var{itmax},@var{tol}])
## PLS1 based on NIPALS
##
##
## Input arguments:
##
## @itemize
## @item
## @code{c}           --- property vector.
## @item
## @code{Y}           --- spectral matrix whose each row is a sample.
## @item
## @code{nPrinComp}   --- number of principal component number.
## @item
## @code{itmax}       --- maximum iteration.
## @item
## @code{tol}         --- convergence criterion.
## @end itemize
## 
## Return values
##
## @itemize
## @item
## @code{W} --- weight matrix.
## @item
## @code{P} --- spectral factors. 
## @item
## @code{T} --- spectral scores.
## @item
## @code{U} --- property scores.
## @item
## @code{b} --- coefficients of inner linear relationships. 
## @item
## @code{q} --- property factors.
## @end itemize
##
## @seealso{functions}
## @end deftypefn

## Author:  Gan, F.
## create date:     2017-10-20
## latest revision: 2017-10-20


function [W,P,T,U,b,q] = plsNipals(c,Y,nPrinComp,itMax,tol)

  if nargin < 5
    error('Please see demo.');
  endif
  if (size(c,2) > 1)
    error('This is for PLS1.');
  endif
  [mRows,nCols] = size(Y);
  nPrinComp = min([mRows,nCols,nPrinComp]);
  U = zeros(mRows,nPrinComp);
  T = ones(mRows,nPrinComp);
  W = zeros(nCols,nPrinComp);
  P = zeros(nPrinComp,nCols);
  q = zeros(nPrinComp,1);
  b = zeros(nPrinComp,1);
  itNum = 0;
  T_old = zeros(mRows,nPrinComp);
  for i = 1:nPrinComp
    U(:,i) = c;
    residualTi = T_old(:,i) - T(:,i);
    while (sum(abs(residualTi)) > tol)
      T_old(:,i) = T(:,i);
      W(:,i) = Y' * U(:,i);
      W(:,i) = W(:,i) / sqrt(W(:,i)' * W(:,i));
      T(:,i) = Y * W(:,i);
      q(i,1) = T(:,i)' * c / (T(:,i)' * T(:,i));
      U(:,i) = c / q(i,1);
      residualTi = T_old(:,i) - T(:,i);
      itNum = itNum + 1;
      if (itNum > itMax)
        printf("Reach maximum iteration at %ith principal component.\n",i);
        break;
      endif
    endwhile
    P(i,:) = T(:,i)' * Y / (T(:,i)' * T(:,i));
    T(:,i) = T(:,i) / sqrt(P(i,:) * P(i,:)');
    W(:,i) = W(:,i) / sqrt(P(i,:) * P(i,:)');
    P(i,:) = P(i,:) / sqrt(P(i,:) * P(i,:)');
    b(i,1) = (U(:,i)' * T(:,i)) / (T(:,i)' * T(:,i));
    Y = Y - T(:,i) * P(i,:);
    c = c - b(i,1) * T(:,i) * q(i,1);
    itnum = 0;
  endfor

endfunction
%!demo
%! load ./Data/CornModelSp.dat;
%! load ./Data/CornModelProp.dat;
%! pcn = 5;
%! itmax = 100;
%! tol = 1e-6;
%! c = CornModelProp - mean(CornModelProp);
%! Y = CornModelSp - repmat(mean(CornModelSp),size(CornModelSp,1),1);
%! [W,P,T,U,b,q] = plsNipals(c,Y,pcn,itmax,tol);
%! figure(1),clf('reset');
%! plot(T(:,1),T(:,2),'o');


