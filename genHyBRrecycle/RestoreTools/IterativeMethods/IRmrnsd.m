function [x, IterInfo] = IRmrnsd(A, b, options)
%
%      [x, IterInfo] = IRmrnsd(A, b, options;
%
%  Modified Residual Norm Steepest Descent
%  Nonnegatively constrained steepest descent method.
%
%   Input: A  -  object defining the coefficient matrix.
%          b  -  Right hand side vector.
% 
%   Optional Intputs:
%     options - Structure that can have:
%                x0      - initial guess; default is x0 = A'*b
%                x_true  - true solution, for testing purposes
%                MaxIter - integer specifying maximum number of iterations;
%                          default is min(size(A), 100)
%                Rtol    - stopping tolerance for the relative residual,
%                                norm(b - A*x)/norm(b)
%                          default is 1e-6
%                NE_Rtol - stopping tolerance for the relative residual,
%                                norm(A'*b - A'*A*x)/norm(A'*b)
%                          default is 1e-6
%                IterBar - 'off' or 'on'
%                          can be used to turn off/on waitbar that shows
%                          iteration progress; default is 'on'
%
% 
%   Output:
%          x  -  solution
%     IterInfo  -  structure containing some information about the iteration
%                Iter     -  actual number of iterations performed
%                Rnrm     -  norm of the residual at each iteration
%                NE_Rnrm  -  norm of the residual at each iteration
%                Xnrm     -  norm of the solution at each iteration
%                Enrm     -  if input options specifies x_true, this contains
%                            the norm of the true relative error at each 
%                            iteration, norm(x - x_true)/norm(x_true)
%                StopFlag -  integer that indicates reason for stopping
%                            iteration:
%                               1 ==> Rtol satisfied
%                               2 ==> NE_Rtol satisfied
%                               3 ==> MaxIter reached
%
% J. Nagy, August, 2011
%
%  References:
%   [1]  J. Nagy, Z. Strakos.
%        "Enforcing nonnegativity in image reconstruction algorithms"
%         in Mathematical Modeling, Estimation, and Imaging, 
%         David C. Wilson, et.al., Eds., 4121 (2000), pg. 182--190. 
%   [2]  L. Kaufman.
%        "Maximum likelihood, least squares and penalized least squares
%         for PET",
%         IEEE Trans. Med. Imag. 12 (1993) pp. 200--214.

% Initialization
defaultopt = struct('x0', 'none', 'MaxIter', [] ,'x_true', 'off', ...
  'Rtol', 1e-6, 'NE_Rtol', 1e-6, 'IterBar', 'on');

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 2
  error('IRmrnsd: Not Enough Inputs')
elseif nargin < 3
  options = [];
end
if isempty(options)
  options = defaultopt;
end

%
% Get options, and make sure everthing is initialized:
%
[m,n] = size(A);
defaultopt.MaxIter = min([m, n, 100]);
options = IRset(defaultopt, options);

x = IRget(options, 'x0', [], 'fast');
MaxIter = IRget(options, 'MaxIter', [], 'fast');
x_true = IRget(options, 'x_true', [], 'fast');
Rtol = IRget(options, 'Rtol', [], 'fast');
NE_Rtol = IRget(options, 'NE_Rtol', [], 'fast');
IterBar = IRget(options, 'IterBar', [], 'fast');

% 
%  Make sure b is shaped as a vector, and not as an image, if the
%  A is full or sparse
%
if isa(A, 'single') || isa(A, 'double')
  [mb, nb] = size(b);
  if size(b,1) ~= m
    b = b(:);
  end
  if size(b,1) ~= m
    error('IRmrnsd: Sizes of input matrix and vector are not compatible')
  end
end

Rtol = norm(b(:)) * Rtol;
trAb = A'*b;
NE_Rtol = norm(trAb(:)) * NE_Rtol;
notrue = strcmp(x_true,{'off'});
if strcmp(x,{'none'}), x = trAb; end

noIterBar = strcmp(IterBar,{'off'});

Rnrm = zeros(MaxIter+1, 1);
Xnrm = zeros(MaxIter+1, 1);
NE_Rnrm = zeros(MaxIter+1, 1);

tau = sqrt(eps);
sigsq = tau;
minx = min(x(:));
%
%  if initial guess has negative values, compensate
%
if minx < 0
     x = x - min(0,minx) + sigsq;
end

if ~notrue
  Enrm = zeros(MaxIter+1, 1);
  nrm_x_true = norm(x_true(:));
  x_error = x(:) - x_true(:);
  Enrm(1) = norm(x_error(:)) / nrm_x_true;
end

r = b - A*x;

Rnrm(1) = norm(r(:));
Xnrm(1) = norm(x(:));

g = -(A'*r);
xg = x .* g;
gamma = g(:)' * xg(:);

NE_Rnrm(1) = sqrt(gamma);

if ~noIterBar
  h = waitbar(0, 'Beginning IRmrnsd iterations: please wait ...');
end

for k = 1:MaxIter
  if Rnrm(k) <= Rtol
    % stop because residual satisfies ||b-A*x||<= Rtol
    StopFlag = 1;
    break
  end
  if NE_Rnrm(k) <= NE_Rtol
    % stop because normal equations residual satisfies ||A'*b-A'*A*x||<= NE_Rtol
    StopFlag = 2;
    break
  end
  if ~noIterBar
    waitbar(k/MaxIter, h)
  end
  
  s = - x .* g;

  u = A*s;
  
  theta = gamma / (u(:)'*u(:));
  neg_ind = s < 0;

  alpha = min( theta, min( -x(neg_ind) ./ s(neg_ind) ) );
  if isempty(alpha)
    alpha = theta;
  end

  x = x + alpha*s;

  z = A'*u;

  g = g + alpha*z;
  xg = x .* g;
  gamma = g(:)' * xg(:);

  if ~notrue
    x_error = x(:) - x_true(:);
    Enrm(k+1) = norm(x_error(:)) / nrm_x_true;
  end

  Rnrm(k+1) = norm(r(:));
  Xnrm(k+1) = norm(x(:));
  NE_Rnrm(k+1) = sqrt(gamma);
end
if k == MaxIter
  % Stop because max number of iterations reached
  StopFlag = 3;
else
  k = k - 1;
end
if ~noIterBar, close(h), end
if isa(A, 'single') || isa(A, 'double')
  x = reshape(x, mb, nb);
end
if nargout==2
  IterInfo.Iter = k;
  IterInfo.Rnrm = Rnrm(1:k+1);
  IterInfo.NE_Rnrm = NE_Rnrm(1:k+1);
  IterInfo.Xnrm = Xnrm(1:k+1);
  if ~notrue
    IterInfo.Enrm = Enrm(1:k+1);
  end
  IterInfo.StopFlag = StopFlag;
end


