## Copyright (C) 2015, Feng Gan <cesgf@mail.sysu.edu.cn; sysucesgf@163.com>
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
## @deftypefn {Function File} {[@var{S}, @var{t}] =} keyspectra (@var{X})
## Generate key spectra and determine the number of component automatically 
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X} is a spectral matrix whose each row is a spectrum.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{S} is a spectral matrix of the key set.
## @item
## @code{t} is the retention time of the key set.
## @end itemize
## @end deftypefn
##
## References:
## - F. Gan, Q.S. Xu, Y. Z. Liang, Two novel procedures for automatic resolution of two-way data from coupled chromatography, Analyst. 126 (2001) 161-168.
## - F.C. Sanchez, J. Toft, B. Van Den Bogaert, D.L. Massart, Orthogonal projection approach applied to peak purity assessment., Anal. Chem. 68 (1996) 79-85. doi:10.1021/ac950496g.
##
## Author: Gan,F.
## Lastest Revision: 2015-07-27
## Create date:      2000-01-16

function [S,t] = keyspectra(X)

  [m, n] = size(X);
  points = ceil(m/2);

  S = [];
  t = [];
  det_i = [];         
  all_det_i = [];     

  init_spectrum = mean(X);
  for j = 1:points
    for i = 1:m
      Y = [init_spectrum; X(i,:)];
      Z = Y*Y';
      det_i = [det_i; det(Z)];
    end
    [maxdetval,maxdetlocat] = max(det_i);           
    all_det_i = [all_det_i det_i/norm(det_i)];
    S = [S; X(maxdetlocat,:)];
    init_spectrum = S;
    det_i = [];
  end

  coef = [];
  for i = 1:points-1
    coef = [coef; all_det_i(:,i)' * all_det_i(:,i+1)];
  end
  difcoef = [];
  for i = 1:points-2
    difcoef = [difcoef (coef(i+1) - coef(i))];
  end
  [maxdifcoef,num] = max(difcoef);
  S = S(1:num,:);
  C = X * pinv(S);
  [maxc, t] = max(C);
  [t, maxt] = sort(t);

  for i = 1:length(t)
    S(i,:) = X(t(i),:) / norm(X(t(i),:));
  end
  S = S';

endfunction

%!demo
%! load ./Data/masartdata.mat
%! [S,t] = keyspectra(X);

