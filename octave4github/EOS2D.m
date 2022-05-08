function [A,B] = EOS2D(X,N,MaxIter,Tol,W)
% Solve X=A*B' by EOS.
%   [A,B] = EOS2D(X,N,MaxIter,Tol,W)
%   or
%   [A,B] = EOS2D(X,N,MaxIter,Tol)
% 
% where
%   - X: Measured data.
%   - N: Number of facto.
%   - MaxIter: Maximum of iterion steps.
%   - Tol: COnvergent criterion. The program convergences 
%          if the relative error of current error to the 
%          error of the last TOL(2) steps is less than Tol(1).
%   - W: Weight of X. Default value is 1.         
% 
% Jihong Wang
% Nov 09, 2000.


if nargin == 4, W = 1; end

SearchSet(1:2:200) = 2.^(-[1:100]);
SearchSet(2:2:201) = -2.^(-[1:100]);
CGRestart = [10,80,1,1,2,1];

nsvdtol = 1e-5;
nsvditermax = 300;    							
[A,S,B] = nsvd(X,N,nsvdtol,nsvditermax);

A = abs(A*sqrt(S));
B = abs(B*sqrt(S));

W = W.*ones(size(X));
At = 0*A; Bt = 0*B;
Ac = ones(size(A)); Bc = ones(size(B)); 
Au = 1./max((W.^2)*B.^2,1); Bu = 1./max((W.^2)'*A.^2,1);
ro = 0;
R = X-A*B';
Q = W(:)'*R(:).^2;

disp(sprintf('\nEOS began at %s\nInitial estimation: Q= %.10e',datestr(now),Q));
disp(sprintf('%23s%.10e ',' ',Q));

QQ = zeros(Tol(2)+6,1); k = Tol(2)+1; 
for i = 1:MaxIter
   
   Ag = (W.*R)*B; Bg = (W.*R)'*A;
   Az = Ac.*Au.*Ag; Bz = Bc.*Bu.*Bg;
   g2 = sum(Az(:).^2)+sum(Bz(:).^2);
   rotmp = Ag(:)'*Az(:)+Bg(:)'*Bz(:);
   beta = 0; 
   if ro ~= 0, beta = rotmp/ro; end; 
   ro = rotmp;
   
   At = beta*At+Az; Bt = beta*Bt+Bz;
   Xv = At*B' + A*Bt';
   ti = W(:)'*Xv(:).^2;
   wi = rotmp;
   alfa = wi/ti;
   Anew = A+alfa*At; Bnew = B+alfa*Bt;
   Anewc = max(Anew,0); Bnewc = max(Bnew,0);
   Rnew = X-Anewc*Bnewc';  
   Qnew = W(:)'*Rnew(:).^2;
   
   Ac = max((Anew<0).*Ac/32,1e-12)+ min((Anew>=0).*Ac*2,1);
   Bc = max((Bnew<0).*Bc/32,1e-12)+ min((Bnew>=0).*Bc*2,1);
   At = max((Anew<0).*At/32,1e-12)+(Anew>=0).*At;
   Bt = max((Bnew<0).*Bt/32,1e-12)+(Bnew>=0).*Bt;
   
   BadStep = 1;
   if Qnew < Q
      BadStep = 0;
      Q = Qnew; 
      R = Rnew;
      A = Anewc; B = Bnewc;
   else
      for j = SearchSet
         Anewc = max(A+j*alfa*At,0);  
         Bnewc = max(B+j*alfa*Bt,0); 
         Rnew = X-Anewc*Bnewc';        
         Qnew = W(:)'*Rnew(:).^2;
         if Qnew < Q
            A = Anewc; B = Bnewc; 
            R = Rnew; 
            Q = Qnew; 
            break; 
         end
      end
   end
   
   disp(sprintf('  %4d  Q= %.10e ||z||= %e ',i,Q,g2));
   QQ(end+1) = Q;
   drmse = abs(mean(QQ([end-Tol(2)]:[end-1])) - Q)/Q;
   if drmse < Tol(1), disp('Stop by deltaQ meets.'); break; end
   
   if [i-k > prod(CGRestart([1,3])) | BadStep == 1] & ...
         [QQ(k)+2*Q-(1+2)*QQ(floor((end-k)/2)) > 0 | ...
            i-k > prod(CGRestart([2,3]))]
      k = i;
      CGRestart(3:6) = CGRestart([4:6,3]);
      Au = 1./max((W.^2)*B.^2,1); Bu = 1./max((W.^2)'*A.^2,1);
      ro = 0;
      Ac = ones(size(Ac)); Bc = ones(size(Bc));
      At = 0*At; Bt = 0*Bt;
   end

#   plot(A)
#   pause
   
end

NB = sqrt(sum(B.^2));
B = B*diag(NB.^(-1));
A = A*diag(NB);

disp(sprintf('\nEOS ended at %s',datestr(now)));
