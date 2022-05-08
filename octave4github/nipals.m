## Copyright (C)  Feng Gan <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
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
## @deftypefn {Function File} {[@var{T},@var{P}] = } nipals (@var{X}, @var{pcn}, @var{itmax}, @var{tol})
## NIPALS algorithm to solve model: X = TP'.
##
## INPUT
## @itemize
## @item
## @code{X}     is a spectral matrix whose each row is a spectrum.
## @item
## @code{pcn}   is the potential number of principal component.
## @item
## @code{itmax} is the maximum iteration number.
## @item
## @code{tol}   is the convergence criterion.
## @end itemize
##
## RETURN
## @itemize
## @item
## @code{T} is the score matrix.
## @item
## @code{P} is the loading matrix.
## @end itemize
##
## @end deftypefn

## Author:  Feng Gan
## Latest revision: 2015-10-01
## Create:          1999-02-21

function [T,P] = nipals(X, pcn, itmax, tol)

  if nargin < 4
    error('Please see demo.');
  endif
  
  [m,n] = size(X);
  pcn = min([m n pcn]); 
  X = X - repmat(mean(X),m,1);
  T = zeros(m,pcn);
  P = zeros(n,pcn);

  for i = 1:pcn

    itnum = 0;
    rho = 1;
    u_j = randn(m,1); 
    p = [];

    while ((rho > tol) && (itnum < itmax))
      p = X' * u_j;
      p = p / norm(p);
      u = X * p;
      rho = u - u_j;
      rho = rho' * rho;
      u_j = u;
      itnum ++;
      if (itnum > itmax)
        disp('Iteration number reaches the maximum value.')
      endif
    endwhile

    T(:,i) = u;
    P(:,i) = p;
    X = X - u * p';

  endfor

endfunction

%!demo
%! load ./Data/smcrm.mat
%! pcn = 3;
%! tol = 1e-5;
%! itmax = 100;
%! [T,P] = nipals(X,pcn,itmax,tol);



