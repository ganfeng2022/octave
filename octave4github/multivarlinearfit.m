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
## @deftypefn {Function File} {[hat_beta,s_e,R2,Ra2] = } multivarlinearfit (c,Y,alpha)
## Multivariate Linear Fitting for c = Yb + e
##
## Input
## @itemize
## @item
##    c     --- response variable
## @item
##    Y     --- predictor variables
## @item
##    alpha --- significant level.
## @end itemize
## Output
## @itemize
## @item
##    hat_b --- estimated beta
## @item
##    s_e   --- standard error of response
## @item
##    R2    --- correlation coefficient
## @item
##    Ra2   --- adjusted R2
## @end itemize
## 
## @end deftypefn

## Author:  Feng Gan
## Latest Revision: 2015-10-12
## Create date: 2015-07-09


function [hat_beta,s_e,R2,Ra2] = multivarlinearfit(c,Y,alpha)

  if (nargin < 2 || nargin > 3)
    error("Input variable must be two or three.");
  endif

  if (!isvector(c))
    error("c must be a numeric vector.");
  endif
  if (!ismatrix(Y))
    error ("Y must be a numberic matrix.");
  endif

  if (rows(c) != rows(Y))
    error("c and Y must contain the same number of rows.");
  endif

  if (nargin < 3)
    alpha = 0.05;
  elseif (!isscalar(alpha))
    error("alpha must be a scalar value")
  endif

  m = rows(Y);
  n = columns(Y) - 1;

  [Q,R] = qr(Y,0);  
  hat_beta = inv(R) * Q' * c;
  hat_c = Y * hat_beta;
  e = c - hat_c;
  s_e = sqrt(sum(e.*e)/(m - n - 1));
  printf("\nThe multivariate linear regression equation is:\n");
  eq_string = "c = ";
  for i = 1:length(hat_beta)
    if i < 2
      eq_string = [eq_string num2str(hat_beta(1))];
    else
      if hat_beta(i) > 0
        eq_string = [eq_string ' + ' num2str(hat_beta(i),"%0.3f") 'y_' num2str(i-1)];
      else
        eq_string = [eq_string ' - ' num2str(abs(hat_beta(i)),"%0.3f") 'y_' num2str(i-1)];
      endif
    endif
  end
  printf(eq_string); printf("\n");
  printf("The standard error (s.e.) = %0.3f\n",s_e);

  var_hat_beta = diag(s_e^2 * inv(Y'*Y));
  se_hat_beta = sqrt(var_hat_beta);
  t_hat_beta = hat_beta ./ se_hat_beta;
  df = m - n - 1;
  p_val = 1 - tcdf(abs(t_hat_beta),df); 
  p_val = 2.0 * p_val;

  R2 = 1 - sum((c - hat_c) .* (c - hat_c))/sum((c - mean(c)) .* (c - mean(c)));
  Ra2 = 1 - (m - 1) * (1 - R2) / (m - n - 1);

  var_string = {'Const.' 'y'};
  printf("\nMultivariate linear regression Output\n");
  printf("----------------------------------------------------------\n");
  printf("Variable    Coeff.    s.e.        t         p-value \n");
  printf("----------------------------------------------------------\n");
  for i = 1:length(hat_beta)
    if i == 1
      printf("%s      %3.3f     %3.3f     %3.3f     %3.3f     %3.3f", var_string{1}, hat_beta(i), se_hat_beta(i), t_hat_beta(i), p_val(i));
    else
      printf("\n%s         %3.3f      %3.3f     %3.3f     %3.3f     %3.3f", [var_string{2} '_' num2str(i-1)], hat_beta(i), se_hat_beta(i), t_hat_beta(i), p_val(i));
    endif
  end
  printf("\n------------------------------------------------------------\n");
  printf("m = %d  R^2 = %2.3f  Ra^2 = %2.3f  s.e. = %2.3f   d.f. = %d \n",m,R2,Ra2,s_e,df);
  printf("------------------------------------------------------------\n\n");

endfunction

%!demo
%! load ./Data/cement.mat;
%! c = Y(:,end);
%! X = [ones(rows(Y),1) Y(:,2:end-1)];
%! [hat_beta,s_e,R2,Ra2] = multivarlinearfit(c,X,0.05);





