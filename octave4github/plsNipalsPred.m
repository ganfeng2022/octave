
## Copyright (C)  <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
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
## @deftypefn {Function File} {[@var{hatc}]} = plsNipalsPred ([@var{X},@var{W},@var{P},@var{b},@var{q},@var{n},@var{barY},@var{barc}])
## Prediction based on the model of PLS1
##
##
## Input arguments:
##
## @itemize
## @item
## @code{X}      --- spectral matrix of unknown samples.
## @item
## @code{W}      --- weight matrix.
## @item
## @code{P}      --- spectral factors.
## @item
## @code{b}      --- coefficients of inner linear relationships.
## @item
## @code{q}      --- property factors.
## @item
## @code{n}      --- principal component number.
## @item
## @code{barY}  --- mean Y, $\bar Y$.
## @item
## @code{barc}  --- mean c, $\bar c$.
## @end itemize
## 
## Return values
##
## @itemize
## @item
## @code{hatc}  --- calculated property,$\hat c$.
## @end itemize
##
## @seealso{functions}
## @end deftypefn

## Author:  Gan, F.
## create date:     2017-10-20
## latest revision: 2017-10-20


function [hatc] = plsNipalsPred(X,W,P,b,q,nPrinComp,barY,barc)

  if nargin < 8
    error('Please see demo');
  endif
  
  [mRows,nCols] = size(X);

  X = X - repmat(barY,mRows,1);
  hatc = zeros(mRows,1);
  T = zeros(mRows,nPrinComp);
  
  for i = 1:nPrinComp
    T(:,i) = X * W(:,i);
    X = X - T(:,i) * P(i,:);
    hatc = hatc + b(i) * T(:,i) * q(i);
  endfor
  
  hatc = hatc + repmat(barc,mRows,1);

endfunction

%!demo
%! load ./Data/CornModelSp.dat;
%! load ./Data/CornModelProp.dat;
%! Y_mc = []; bar_Y = []; c_mc = []; bar_c = [];
%! [Y_mc,bar_Y] = meancent(CornModelSp,'col'); 
%! [c_mc,bar_c] = meancent(CornModelProp,'col');
%! pcn = 5;
%! itmax = 100;
%! tol = 1e-6;
%! [W,P,T,U,b,q] = plsNipals(c_mc,Y_mc,pcn,itmax,tol);
%! [hat_c] = plsNipalsPred(CornModelSp,W,P,b,q,pcn,bar_Y,bar_c)
%! figure(1),clf('reset');
%! plot(CornModelProp,hat_c,'o');
%! xlabel('Real');ylabel('Calculated')


