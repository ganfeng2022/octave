
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
## @deftypefn {Function File} {[@var{forwardEigenValues}, @var{reverseEigenValues}] =} efa (@var{X}, @var{nEigenValues})
## Evolving factor analysis
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X}                  --- a spectral matrix whose each row is a spectrum.
## @item
## @code{nEigenValues}       --- the number of eigen values to be collected.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{forwardEigenValues} --- the forward eigen values.
## @item
## @code{reverseEigenValues} --- the reverse eigen values.
## @end itemize
## @end deftypefn
##
## References:
## - Gampp, H., Maeder, M., Meyer, C. J. & Zuberbuhler, A. D. Talanta 32, 95–101 (1985).
##
## Author:  Feng Gan
## Lastest Revision: 2017-10-28
## Create date:      2008-07-08


function [forwardEigenValues,reverseEigenValues] = efa(X, nEigenValues);

  if (nargin < 2)
    nEigenValues = 5;  
  end

  [mRows, nCols] = size(X);
  if nEigenValues > (min([mRows nCols]))/2
    error('You set a too big nev value!');
  end

  forwardEigenValues = zeros(mRows - 1, nEigenValues);
  reverseEigenValues = zeros(mRows - 1, nEigenValues);

  for i = 2:mRows
    subX = X(1:i,:);
    mRowsSub = size(subX,1);
    if mRowsSub < nCols
      subX = subX' * subX;   
    else
      subX = subX * subX';
    endif
    eigenValues = svd(subX);
    eigenValues = sqrt(eigenValues(1:min(mRowsSub, nEigenValues)));
    forwardEigenValues(i-1, 1:min(mRowsSub, nEigenValues)) = eigenValues'; 
  end

  for i = (mRows-1):-1:1
    subX = X(i:mRows, :);
    mRowsSub = size(subX,1);
    if mRowsSub < nCols
      subX = subX' * subX;
    else
      subX = subX * subX';
    end
    eigenValues = svd(subX);
    eigenValues = eigenValues(1:min(mRowsSub,nEigenValues));
    eigenValues = sqrt(eigenValues(1:min(mRowsSub,nEigenValues)));
    reverseEigenValues(i, 1:min(mRowsSub,nEigenValues)) = eigenValues';
  end

  for i = 1:mRows-1
    for j = 1:nEigenValues
      if forwardEigenValues(i,j) > 0
        forwardEigenValues(i,j) = log10(forwardEigenValues(i,j));
      end
      if reverseEigenValues(i,j) > 0
        reverseEigenValues(i,j) = log10(reverseEigenValues(i,j));
      end
    end
  end

  for i = 1:nEigenValues
    forwardEigenValues(1:i,i) = forwardEigenValues(i,i);   
    reverseEigenValues(mRows-i:mRows-1,i) = reverseEigenValues(mRows-i,i); 
  end

endfunction

%!demo
%! load ./Data/masartdata.mat;
%! [fev,rev] = efa(X, 6);
%! plot(fev)
%! hold on
%! plot(rev,'-.')



