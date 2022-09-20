function [X,rho,eta,F] = cglsAIR(A,b,K)
%CGLSAIR Conjugate gradient algorithm applied implicitly to the normal equations
%
% X = cgls(A,b,K)
%
% Performs max(K) steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
%
% The function returns the iterates associated with the values in the
% vector K, stored as columns of the matrix X.  The function is a simplified
% version of the function cgls from Regularization Tools.

% References: A. Bjorck, "Numerical Methods for Least Squares Problems",
% SIAM, Philadelphia, 1996.

% Per Christian Hansen, IMM, May 22, 2012.

% Initialization.
k = max(K);
if (k < 1), error('Number of steps k must be positive'), end
n = size(A,2); X = zeros(n,length(K));

% Prepare for CG iteration.
x = zeros(n,1);
d = A'*b;
r = b;
normr2 = d'*d;

% Iterate.
ksave = 0;
for j=1:k

  % Update x and r vectors.
  Ad = A*d; alpha = normr2/(Ad'*Ad);
  x  = x + alpha*d;
  r  = r - alpha*Ad;
  s  = A'*r;

  % Update d vector.
  normr2_new = s'*s;
  beta = normr2_new/normr2;
  normr2 = normr2_new;
  d = s + beta*d;
  
  % Save, if wanted.
  if any(K==j)
      ksave = ksave + 1;
      X(:,ksave) = x;
  end

end
