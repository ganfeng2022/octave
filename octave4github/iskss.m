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
## @deftypefn {Function File} {[@var{C}, @var{S}, @var{Res}] =} iskss (@var{X}, @var{siglev})
## Improved Stepwise Key Spectra Selections 
## with the model @code{X = CS}.
##
## @itemize
## @item
## @code{X} is a spectral matrix whose each row is a spectrum.
## @item
## @code{siglev} is a significant level for eliminating negative part.
## @end itemize
##
## Return values
##
## @itemize
## @item
## @code{C} is the concentration matrix
## @item
## @code{S} is the spectral matrix of the key set.
## @item
## @code{Res} is the residual matrix
## @end itemize
## @end deftypefn
##
## References:
## - F. Gan, Q.S. Xu, Y. Z. Liang*, Two novel procedures for automatic resolution of two-way data from coupled chromatography, Analyst. 126 (2001) 161-168.
## - L.M. Peng, F. Gan*, H.H. Kong, X.D. Lin, and C.L. Chen, Resolution of Overlapping Peaks from Hyphenated Chromatography by a Procedure Combining Improved Multivariate Curve Resolution Techniques (manuscript)
##
## Author: Gan,F.
## Lastest Revision: 2015-07-28
## Create date:      2011-08-26

function [C,S,Res] = iskss(X,siglev);

  if (nargin < 2)
    siglev = 3.0; ## 缺省为3倍背景的标准偏差。需根据具体情况进行调整。
  endif
  X_sub = [];		
  s_i_ks = [];    
  C_tmp = [];		
  S_tmp = [];		
  U = [];
  s = [];		    
  V = [];
  c_val = 0.0; 			
  c_bk = [];			    
  locat = 0;	            
  V_ks_right = [];		
  V_ks_left = [];
  mx = 0;
  nx = 0;
  pcnum = 0;
  [S,t] = keyspectra(X);
  [mx,nx] = size(X);
  pcnum = length(t);
  if pcnum < 3
    error('You do not need me to resolute the data.');
  end
  c = zeros(mx,pcnum);					
  s_i_ks = X(t(1),:);
  s_i_ks = s_i_ks/norm(s_i_ks);
  for i = t(2)+1:mx
    X_sub = X(i:mx,:); 					
    [U,s,V] = svd(X_sub);				
    S_tmp = [s_i_ks; V(:,1:pcnum-1)'];			
    C_tmp = X * pinv(S_tmp);					
    c_bk = C_tmp(i:mx, 1);				
    [c_val,loca] = min(c_bk);			
    if (abs(c_val) < siglev*std(c_bk))			
      c(:,1) = C_tmp(:,1);			
      break;
    end 
  end

  s_i_ks = X(t(pcnum),:);
  s_i_ks = s_i_ks/norm(s_i_ks);
  for i = t(pcnum-1)-1:-1:1
    X_sub = X(1:i,:);
    [U,s,V] = svd(X_sub);
    S_tmp = [V(:,1:pcnum-1)';s_i_ks];			
    C_tmp= X*pinv(S_tmp);
    c_bk = C_tmp(1:i, pcnum);
    [c_val,loca] = min(c_bk);
    if (abs(c_val) < siglev*std(c_bk))
      c(:,pcnum) = C_tmp(:,pcnum);
      break;
    end
  end 

  for i = 2:pcnum-1
    s_i_ks = X(t(i),:);						
    s_i_ks = s_i_ks/norm(s_i_ks);			
    X_sub = X(t(i+1):mx,:);					
    [U,s,V] = svd(X_sub);
    V_ks_right = V(:,1:pcnum-i);			
    locat = 0;
    for j = t(i-1):-1:1
      X_sub = X(1:j,:);							
      [U,s,V] = svd(X_sub);
      S_tmp=[V(:,1:i-1)';V_ks_right';s_i_ks];		
      C_tmp = X*pinv(S_tmp);
      c_bk = C_tmp(1:j-1,pcnum);
      [c_val,loca] = min(c_bk);
      if (abs(c_val) < siglev*std(c_bk))
        locat = j;
        break;
      end
    end
    X_sub = X(1:locat,:);		
    [U,s,V] = svd(X_sub);
    V_ks_left = V(:,1:i-1);		
    for j = t(i+1):mx
      X_sub = X(j:mx,:);
      [U,s,V] = svd(X_sub);
      S_tmp = [V_ks_left';V(:,1:pcnum-i)';s_i_ks];
      C_tmp = X*pinv(S_tmp);
      c_bk = C_tmp(j+1:mx,pcnum);
      [c_val,loca] = min(c_bk);
      if (abs(c_val) < siglev*std(c_bk))
        c(:,i) = C_tmp(:,pcnum);
        break;
      end
    end
  end
  C = c;
  S = pinv(C)*X;
  Res = X - C*S;
  S = S';
endfunction

%!demo
%! load ./Data/masartdata.mat
%! siglev = 2.5;
%! [C,S,Res] = iskss(X,siglev);
%! figure(1),clf('reset');
%! plot(C);
%! figure(2),clf('reset');
%! plot(S);
%! figure(3),clf('reset');
%! plot(Res);
