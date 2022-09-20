function [Q,k,Qr] = getcovfun(nvec,type, nu, ell, theta)
%
% Get a spatial covariance matrix
%
% nvec is [n n] or [n,n,nt]

xmin = zeros(size(nvec));           % Coordinates of left corner
xmax = ones(size(nvec));           % Coordinates of right corner
if nargin < 5
    theta = ones(size(nvec));
end
sigma2 = 1;
switch type
    case 'matern'
        k = @(r)sigma2*matern(r,nu,ell);
    case 'rational'
        k = @(r) rational_quadratic(r,nu,ell);
    case 'exponential'
        k = @(r) exp(-((r/nu).^2)/2);
end

% Define prior covariance matrix
Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x) toeplitzproduct(x, Qr, nvec);
if length(nvec) == 2
    Q = funMat(Qfun,Qfun,[nvec(1)*nvec(2) nvec(1)*nvec(2)]);
elseif length(nvec) == 3
    Q = funMat(Qfun,Qfun, [nvec(1)*nvec(2)*nvec(3) nvec(1)*nvec(2)*nvec(3)]);
end
