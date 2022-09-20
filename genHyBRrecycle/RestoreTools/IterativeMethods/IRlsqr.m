function [x, IterInfo] = IRlsqr(A, b, options)
%
%      [x, IterInfo] = IRlsqr(A, b, options);
%
%  LSQR (conjugate gradient) method for Iterative Image Restoration
%  Note: This is a wrapper function for MATLAB's built-in LSQR function.
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
%
% 
%   Output:
%          x  -  solution
%   IterInfo  -  structure containing some information about the iteration
%                Iter     -  actual number of iterations performed
%                Rnrm     -  norm of the residual at each iteration
%                NE_Rnrm  -  scaled norm of the normal equations residual 
%                            at each iteration
%                Xnrm     -  norm of the solution at final iteration
%                Enrm     -  if input options specifies x_true, this contains
%                            the norm of the true relative error at final 
%                            iteration, norm(x - x_true)/norm(x_true)
%                StopFlag -  integer that indicates reason for stopping
%                            iteration:
%                               0 ==> Rtol satisfied
%                               1 ==> MaxIter reached
%
%                               
%
% J. Nagy, August, 2011
%
%  References: [1] C. Paige and M. Saunders. "LSQR: An algorithm for sparse 
%                  linear equations and sparse least squares",
%                  ACM Trans. Math. Soft.", 8 (1982), pp. 43-71.
%              [2] C. Paige and M. Saunders. "Algorithm 583 LSQR: An 
%                  algorithm for sparse linear equations and sparse least 
%                  squares", ACM Trans. Math. Soft.", 8 (1982), pp. 43-71.

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
  error('IRlsqr: Not Enough Inputs')
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

[mb, nb] = size(b);
% 
%  For LSQR we'll make sure b is a vector
if isa(A, 'single') || isa(A, 'double')
  if size(b(:),1) ~= m
    error('IRlsqr: Sizes of input matrix and vector are not compatible')
  else
      b = b(:);
  end
   
end

notrue = strcmp(x_true,{'off'});
if strcmp(x,{'none'}), x = A'*b; end

noIterBar = strcmp(IterBar,{'off'});
if isa(A, 'single') || isa(A, 'double')
  [x,StopFlag,relres,iter,Rnrm,NE_Rnrm,Enrm] = modified_lsqr(A,b,Rtol,MaxIter,[],[],x(:),x_true,noIterBar);
else
  b = b(:);
  [x,StopFlag,relres,iter,Rnrm,NE_Rnrm, Enrm] = modified_lsqr(@(v, tflag)Amult(v, A, tflag),b,Rtol,MaxIter,[],[],x(:), x_true, noIterBar);
end

Xnrm = norm(x(:));
x = reshape(x, mb, nb);
if nargout==2
  IterInfo.Iter = iter;
  IterInfo.Rnrm = Rnrm;
  IterInfo.NE_Rnrm = NE_Rnrm;
  IterInfo.Xnrm = Xnrm;
  if ~notrue
    IterInfo.Enrm = Enrm;
  end
  IterInfo.StopFlag = StopFlag;
end

function Av = Amult(v, A, transp_flag)
%
%  This is needed for a user defined object A for which ctranspose and
%  mtimes have been overloaded.  In particular, for IR codes we need
%  this when using a psfMatrix object.
%

% Apply A or A^T to vector v
sqrtSize_v = sqrt(size(v,1));
v = reshape(v,sqrtSize_v,sqrtSize_v);

if strcmp(transp_flag, 'notransp')==1
    Av = A*v;
elseif strcmp(transp_flag, 'transp')==1
    Av = A'*v;
else
    error('invalid transp_flag option')
end

Av = Av(:);


