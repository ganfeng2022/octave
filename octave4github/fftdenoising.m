
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
## @deftypefn {Function File} {[@var{denoisedSignal}]} = fftdenoising ([@var{noisedSignal}, @var{nSelectedFrequency}])
## denosing using FFT
##
## Input varibles
##
## @itemize
## @item
## @code{noisedSignal}       --- signal with noise
## @item
## @code{nSelectedFrequency} --- number of selected frequency
## @end itemize
## 
## Return values
##
## @itemize
## @item
## @code{denoisedSignal}     --- denoised signal
## @end itemize
##
## @seealso{fft,ifft}
## @end deftypefn

## Author:  Feng GAN
## create date:     2003-07-27
## latest revision: 2016-11-04

function [denoisedSignal] = fftdenoising(noisedSignal,nSelectedFrequency)

  if (nargin < 2)
    error('Please see demo.');
  endif
  
  frequency = fft(noisedSignal);
  mRows = length(frequency);
  selectedFrequencyRange = ones( mRows, 1 );
  selectedFrequencyRange( nSelectedFrequency : (mRows - nSelectedFrequency - 1) ) = 0;
  retainedFrequency = frequency .* selectedFrequencyRange;
  denoisedSignal = ifft(retainedFrequency);
  denoisedSignal = real(denoisedSignal);

endfunction

%!demo
%! load ./Data/noisedata1.txt
%! x = noisedata1(:,2);
%! sf = 12;
%! [x_new] = fftdenoising(x,sf);
%! plot(x_new)

