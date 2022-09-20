function x_true = dynamic_test(T, n)
[x,y] = meshgrid(-n/2:n/2-1,-n/2:n/2-1);
A = zeros(n);
A(x.^2 + y.^2 < 1000)=1;
%[x,y] = meshgrid(linspace(-1,1,n),linspace(-1,1,n));
k = 1;
x_true = [];
rng(0)
c = linspace(1,50,T);
for k = 1:T
    %       m1 = sin(5*(10^-5));
    %       m2 = sin(7*(10^-5));
    %
    %       sx = 1 + nthroot(5*m1,4).*x + (nthroot(5*m1,4).*x).^2 + (nthroot(5*m1,4).*x).^3 + (nthroot(5*m1,4).*x).^4;
    %       sy = 1 + nthroot(5*m2,4).*y + (nthroot(5*m2,4).*y).^2 + (nthroot(5*m2,4).*y).^3 + (nthroot(5*m2,4).*y).^4;
    %
    %       ssx = sx.*x;
    %       ssy = sy.*y;
    %      x_true(:,:,k) = interp2(x,y,A,ssx,ssy,'spline',0);
    A((x-1).^2 +  (y+1).^2 < c(k)*6)=.5;
    x_true(:,:,k) = A;
end
end