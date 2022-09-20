function [f] = l1_TV(A,b,lambda_1,lambda_TV,N, obj)
% joint l_1 and TV regularization: find the silution to the problem:
%          x=argmin_x>=0
%          1/2||Ax-b||^2+lambda_1||x||_1+lambda_TV*sum_k|x_m(k)-x_n(k)|
%   By substituting the original object function by a separable paraboidal
%   surrogate function that satisfies a set of majorization
%   conditions(majorization) and then minimize the surrogate function(minimization) 
%   we can solve the original minimization problem.
%
% Remark1: the l_1 penalty term here is to promote sparsity and the TV
%        penalty term is to introduce smoothness so as to preserve the edges of an
%        image. The larger the corresponding penalty weight, the stronger the effect
%        of sparsifying or smoothing respectively. 
% Remark2: nonnegativity constraint is also considered in this problem and 
%        nonnegativity is achieved by projecting the reconstruction result
%        onto the nonnegative orthant.
% Remark3: the algorithm starts with a randomly picked positive guess x_0 if no 
%        a priori information of the distribution of the fluorescence
%        intensity is given in advance.
% Remark4: the image vector x is first reshaped 
M=size(A,2);
x=rand(M,1);
ATA=A'*A;
ATb=A'*b;
inv_D=zeros(M,M);
for n=1:N
    
    X=reshape(x,obj.recon_grd);
    C = computeC(X);
    abs_C=abs(C);
    v=C*x;
    z=zeta(v);
    H=2*lambda_TV*abs_C'*z;
    for j=1:M   
        inv_D(j,j)=1/(norm(ATA(:,j),1)+H(j));
    end
    G_phi=ATA*x-ATb+lambda_1*ones(M,1)+lambda_TV*C'*(v.*z);
    xx=max(x-inv_D*G_phi,0);
    x=xx;
end
f=x;
end


function C = computeC(b)


    dim = size(b);
    c = zeros(dim + 2);
    c(2:end-1, 2:end-1, 2:end-1) = b;

    right = 2 * (c(2:end-1,2:end-1,2:end-1) > c(2:end-1,3:end,2:end-1)) -1;
    right(:,end,:) = 0;


    left = 2*(c(2:end-1,2:end-1,2:end-1) >= c(2:end-1,1:end-2,2:end-1)) -1;
    left(:,1,:) = 0;

    top = 2*(c(2:end-1,2:end-1,2:end-1) > c(1:end-2,2:end-1,2:end-1)) -1;
    top(1,:,:) = 0;

    bottom = 2*(c(2:end-1,2:end-1,2:end-1) >= c(3:end,2:end-1,2:end-1)) -1;
    bottom(end,:,:) = 0;

    front = 2*(c(2:end-1,2:end-1,2:end-1) > c(2:end-1,2:end-1,1:end-2)) -1;
    front(:,:,1) = 0;

    back = 2*(c(2:end-1,2:end-1,2:end-1) >= c(2:end-1,2:end-1,3:end)) -1;
    back(:,:,end) = 0;


    hor = left + right;
    ver = top + bottom;
    orth = front + back;

    C = [hor(:)'; ver(:)'; orth(:)'];

end

