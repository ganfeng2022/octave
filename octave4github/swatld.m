
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
## @deftypefn {Function File} {[@var{A},@var{B},@var{C},@var{res}]} = swatld ([@var{X},@var{pnc},@var{K},@var{tol},@{itmax}])
## Self-weighted Alternating Trilinear Decomposition (SWATLD)
##
## Reference:
## Chen, Z.P.; Wu, H. L.; Jiang, J. H.; Li, Y. & Yu, R. Q. A novel trilinear decomposition algorithm for second-order linear calibration Chemometrics and Intelligent Laboratory Systems, 2000 , 52 , 75-86.
##
## Input arguments:
##
## @itemize
## @item
## @code{X}   --- column-wise collapsed three-way array.
## @item
## @code{pnc} --- potential number of component.
## @item
## @code{K}   --- K is the number of measured two-way data array.
## @item
## @code{tol}   --- tolerence of convergence
## @item
## @code{itmax} --- maximum iteration number
## @end itemize
## 
## Return values
##
## @itemize
## @code{A}   --- loading matrix of mode A.
## @item
## @code{B}   --- loading matrix of mode B.
## @item
## @code{C}   --- loading matrix of mode C.
## @item
## @code{res} --- residuals in the iterations.
## @end itemize
##
## @seealso{functions}
## @end deftypefn

## Author:  Gan, F.
## create date:     2018-01-11
## latest revision: 2018-01-11
##
## This program is revised from the swatld.m provided by Prof. Chen,Z.P.
## 
function [A,B,C,res] = swatld(X,pnc,K,tol,itmax)

  ## dimensions
  [m,n] = size(X)
  I = m/K;
  J = n;

  ## initialization
  A = rand(I,pnc);
  B = rand(J,pnc);
  C = rand(K,pnc);
  res = 1.0;
  it = 1;
  err = 1.0
  e = 1e-6;
  while (err > tol && it < itmax)
    tol_0 = res(it);
    ## iteration on K slabs 
    for k = 1:K
      X_k = X(((k-1)*I+1):(k*I),:);
      C_1 = diag(pinv(A,e) * X_k * B)  ./ diag(B'*B);
      C_2 = diag(pinv(B,e) * X_k' * A) ./ diag(A'*A);
      C(k,:) = (C_1' + C_2') / 2.0;
      C(k,:) = max(C(k,:),0);
    endfor
    ## iteration on J slabs
    for j = 1:J
      X_j = reshape(X(:,j),I,K);
      B_1 = diag(pinv(A,e) * X_j * C)  ./ diag(C' * C);
      B_2 = diag(pinv(C,e) * X_j' * A) ./ diag(A' * A);
      B(j,:) = (B_1' + B_2') / 2.0;
      B(j,:) = max(B(j,:));
    endfor
    ## normalization
    B = B * diag(1./(sqrt(diag(B'*B) + (diag(B'*B)==0))));
    ## iteration on I slabs
    for i = 1:I
      X_i = X(i:I:I*K,:);
      A_1 = diag(pinv(C,e) * X_i * B)  ./ diag(B'*B);
      A_2 = diag(pinv(B,e) * X_i' * C) ./ diag(C'*C);
      A(i,:) = (A_1' + A_2') / 2.0;
      A(i,:) = max(A(i,:),0);
    endfor
    ## normalization
    A = A * diag(1./(sqrt(diag(A'*A) + (diag(A'*A)==0))));
    ## calculate the residuals
    tmp = 0;
    for k = 1:K
      tmp1 = A * diag(C(k,:)) * B' - X((k-1)*I+1:k*I,:);
      tmp2 = tmp1 .* tmp1;
      tmp = tmp + sum(tmp2(:));
    endfor
    res = [res tmp];
    it = it + 1;
    err = abs((res(it) - res(it-1)) / res(it-1));

  endwhile

endfunction

%!demo
%! [X_IJK,XIJK_rowwise,XIJK_colwise,A0,B0,C0] = simu3d(80,80,6,3,'fluo3');
%! [A,B,C,res]=swatld(XIJK_colwise,4,18,1e-4,200);
%! figure(1),clf('reset');
%! plot(A)
%! figure(1),clf('reset');
%! plot(A); 
%! figure(2),clf('reset');
%! plot(B)
%! figure(3),clf('reset');
%! plot(C)
%! figure(4),clf('reset');
%! plot(res(2:end))


