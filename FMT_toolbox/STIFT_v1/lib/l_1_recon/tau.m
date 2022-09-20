function [y] = tau(alpha,x)
%shrinkage operator tau_alpha is a transform applied to the given vector x,
%                   it is the main ingredient of the
%                   shrinkage/soft-thresholding step of the class of
%                   iterative shrinkage-thresholding algorithms(ISTA).
%   here, alpha:=lambda*t, where
%              lambda: regularization weight in the minimization problem 
%              t: appropriate stepsize        
d=size(x,1);
y=zeros(d,1);
for i=1:d
    y(i)=max(abs(x(i))-alpha,0)*sign(x(i));
end
end

