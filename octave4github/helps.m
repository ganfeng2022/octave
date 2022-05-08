## Copyright (C) 2015, Feng Gan <cesgf@mail.sysu.edu.cn;sysucesgf@163.com>
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

## -*- teselectXnfo -*-
## @deftypefn {Function File} {[@var{C}, @var{S}, @var{Residuals}] =} helps (@var{X}, @var{nPrinComp}, @var{selectInfor})
## Heuristic evolving latent projections
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X}           --- a two way array whose each row is a spectrum.
## @item
## @code{nPrinComp}   --- the number of components.
## @item
## @code{selectInfor} --- a matrix containing selective regions and zero concentration regions.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{C}         --- the concentration matrix.
## @item
## @code{S}         --- the spectral matrix.
## @item
## @code{Residuals} --- the residual matrix.
## @end itemize
## @end deftypefn
##
## References:
##
#[1] O.M. Kvalheim, Y.Z. Liang, Heuristic evolving latent projections: resolving two-way multicomponent data. 1. Selectivity, latent-projective graph,datascope, local rank, and unique resolution, Anal. Chem. 64 (1992) 936
#
#[2] Y.Z. Liang, O.M. Kvalheim, H.R. Keller, D.L. Massart, P. Kiechle, F. Erni, et al., Heuristic evolving latent projections: resolving two-way multicomponent data. 2. Detection and resolution of minor constituents, Anal. Chem. 64 (1992) 946
##
## Author:  Feng Gan
## Latest revision: 2015-10-17
## Create:          1999-01-20

function [C, S, Residuals] = helps(X, nPrinComp, selectInfor)

  if (nargin < 3)
    error('Please see demo.');
  endif

  C = [];
  S = [];
  Residuals = [];

  mRowsSelectComp = size(selectInfor,1);

  for i = 1:mRowsSelectComp

    [U,L,V] = svd(X);
    U = U(:,1:nPrinComp) * L(1:nPrinComp,1:nPrinComp);

    range1 = selectInfor(i,1):selectInfor(i,2);
    range2 = selectInfor(i,3):selectInfor(i,4);
    range3 = selectInfor(i,5):selectInfor(i,6);

    selectX = X(range2,:);
    [selectU,selectS,selectV] = svd(selectX);
    selectU = selectU * selectS;
    Si = selectV(:,1);
    Ci = selectU(:,1);

    [v,locat] = max(abs(Ci));
    if (Ci(locat) < 0)
      Ci = -1.0 * Ci;
      Si = -1.0 * Si;
    end

    Ci = [zeros(size(range1'));Ci;zeros(size(range3'))];
    Ui = [U(range1',:);U(range2',:);U(range3',:)];
    Ri = inv(Ui' * Ui) * Ui' * Ci;
    Ci = U * Ri;
    C(:,i) = Ci;
    S(:,i) = Si;
    X = X - Ci * Si';
    nPrinComp = nPrinComp - 1;
    
  endfor

  Residuals = X;

endfunction

%!demo
%! load ./Data/masartdata.mat
%! nPrinComp = 5;
%! selectInfor = [1, 10, 11, 16, 22, 50;
%!                1, 30, 36, 43, 44, 50;
%!                1, 23, 32, 35, 36, 50;
%!                1, 19, 29, 32, 33, 50];
%! [C, S, Res] = helps(X, nPrinComp, selectInfor);
%! figure(1),clf('reset')
%! plot(C)
%! figure(2),clf('reset')
%! plot(S)
%! figure(3),clf('reset')
%! plot(Res)


