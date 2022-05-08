
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
## @deftypefn {Function File} {[@var{A},@var{B},@var{C},@var{Res}]} = parafac ([@var{X_IJK},@var{pcn},@var{maxit},@var{tol}])
## Parallel factor analysis (PARAFAC) 
##
## This program is the simplest version of PARAFAC for beginner. The mathematical
## model is X_k = A * C_k * B'. Only nonnegative constraints on the
## all Modes are used.
##
## Input varibles
## @itemize
## @item
## @code{X_IJK}  three-way data array. I and J are measuring chenel
## and K is the number of occations.
## @item
## @code{pcn}    number of response species.
## @item
## @code{maxit}  maximum iteration number.
## @item
## @code{tol}    tolerence of convergence.
## @end itemize
## 
## Return values
## @itemize
## @item
## @code{A}   loading matrix of A, it can be two-way excitation
## fluorence matrix.
## @item
## @code{B}   loading matrix of B, it can be two-way emission
## fluorence matrix.
## @item
## @code{C}   loading matrix of C, it is a concentration matrix. In
## this program, the columns of C is the concentrations of the
## species.
## @item
## @code{Res} residuals
## @end itemize
##
## @seealso{functions
## @end deftypefn

## Author: Gan, F.
## create date:     2017-07-22
## latest revision: 2017-08-02


function [A,B,C,Res] = parafac(X_IJK,pcn,maxit,tol)

  if (nargin < 4)
    error ('Please see demo.');
  endif

  [I,J,K] = size(X_IJK);
  X1 = [];
  X2 = [];
  A = [];
  B = [];
  C = [];
  Res = [];  
  init_tol = sum(X_IJK(:));

  ## string out each slabs of X_IJK horizontally and vectically, respectively.
  for i = 1:K
    X1 = [X1 X_IJK(:,:,i)];
    X2 = [X2;X_IJK(:,:,i)];
  endfor

  ## initialize B and C
  B = rand(J,pcn);
  C = rand(pcn,K);   ## each column is the diagonal of a C_k.
  ## you can also use following initialization for C.
  for k = 1:K
    [u,s,v] = svd(X_IJK(:,:,k));
    ev = diag(s);
    C(:,k) = ev(1:pcn);
    C(:,k) = C(:,k) ./ norm(C(:,k));
  endfor
  
  ## iterative least squares
  for i = 1:maxit

    ## string out each of C_k * B' horizongtally.
    CxB = [];
    for k = 1:K
      CxB = [CxB (diag(C(:,k)) * B')];
    endfor
    ## calculate A
    A = X1 * CxB' * inv(CxB * CxB');
      
    ## string out each of A * C_k vertically.
    AxC = [];
    for k = 1:K
      AxC = [AxC;(A * diag(C(:,k)))];
    endfor
    ## calculate B
    B = inv(AxC' * AxC) * AxC' * X2; 
    B = B';
    
    ## normalize the columns of A and B into unit length.
    for j = 1:pcn
      A(:,j) = A(:,j) ./ norm(A(:,j));
      B(:,j) = B(:,j) ./ norm(B(:,j));
    endfor
    ## calculate C_k
    for k = 1:K
      C_k = pinv(A,tol) * X_IJK(:,:,k) * pinv(B',tol);
      C(:,k) = diag(C_k);
    endfor

    ## nonnagetive constraints.
    A = abs(A); 
    B = abs(B); 
    C = abs(C); 

    ## calculate residuals
    tmp = 0.0;
    for k = 1:K
      delta_X_k = A * diag(C(:,k)) * B' - X_IJK(:,:,k);
      XxX = delta_X_k' * delta_X_k;
      tmp = tmp + sum(XxX(:));      
    endfor
    Res = [Res tmp];
    printf('The residual = %f at %dth iteration.\n', Res(i), i);
    tmp = (Res(i) - init_tol) / init_tol;

    ## evalulate convergence
    if (abs(tmp) < tol)
      break;
    else
      init_tol = Res(i);       
    endif
    
  endfor

endfunction

%!demo
%! load ./Data/X_IJK_1.mat; ## A simulated three-way data array with K slabs.
%! pcn = 3; maxit = 1000; tol = 1e-4; 
%! [A,B,C,Res] = parafac(X_IJK,pcn,maxit,tol);
%! figure(1),clf('reset'),plot(A),title('Loading matrix of mode A');
%! figure(2),clf('reset'),plot(B),title('Loading matrix of mode B');
%! figure(3),clf('reset'),plot(C'),title('Loading matrix of mode C');
%! figure(4),clf('reset'),plot(Res),title('Residuals'),xlabel('Iteration number');
