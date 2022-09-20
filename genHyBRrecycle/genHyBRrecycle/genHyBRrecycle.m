function [x_out, output, trunc_mats] = genHyBRrecycle(A, b, P, Q, R, options, trunc_options, trunc_mats)
%
% [x_out, output, trunc_mats] = genHyBRrecycle(A, b, P, Q, R, options, trunc_options, trunc_mats)
%
% genHyBRrecycle is a Golub-Kahan-based generalized hybrid projection methods that can
% exploit compression and recycling techniques in order to solve a broad
% class of inverse problems where memory requirements or high computational
% cost may otherwise be prohibitive.
%
% Inputs:
%                A : either (a) a full or sparse matrix
%                           (b) a matrix object that performs matrix*vector
%                                 and matrix'*vector operations
%                b : rhs vector
%                P : left preconditioner, P_left, OR
%                  : cell containing left and right preconditioner (optional)
%                          {P_left, P_right}
%                   Note: Preconditioning is not yet implemented for the
%                   recycling process
%             Q, R : covariance matrices
%
% options : structure with the following fields (optional)
%         InSolv - solver for the inner problem: [none | TSVD | {Tikhonov}]
%         RegPar - a value or method to find the regularization parameter:
%                       [non-negative scalar | GCV | {WGCV} | optimal]
%                   Note: 'optimal' requires x_true
%          Omega - if RegPar is 'WGCV', then omega must be
%                       [non-negative scalar | {adapt}]
%           Iter - maximum number of Lanczos iterations:
%                       [ positive integer | {min(m,n,100)} ]
%         x_true - True solution : [ array | {off} ]
%                Returns error norms with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%         BegReg - Begin regularization after this iteration:
%                   [ positive integer | {2} ]
%         Reorth - it is the option for HyBR rather than HyBR recycle
%             Vx - extra space needed for finding optimal reg. parameters
%        FlatTol - Tolerance for detecting flatness in the GCV curve as a
%                    stopping criteria
%                   [ non-negative scalar | {10^-6}]
%         MinTol - Window of iterations for detecting a minimum of the GCV curve
%                    as a stopping criteria
%                   [ positive integer | {3}]
%         ResTol - Residual tolerance for stopping the LBD iterations,
%                    similar to the stopping criteria from [1]: [atol, btol]
%                   [non-negative scalar  | {[10^-6, 10^-6]}]
%       Note: options is a structure created using the function 'HyBRset'
%               (see 'HyBRset' for more details)
%
% trunc_options : structure with additional parameters for truncation
%        nOuter - number of outer iterations
%        nInner - maximum storage of solution vector space
%        max_mm  - maminimum number of vectors to save at compression
%        compress  - method for compression
%
% trunc_mats : recycled matrices, W, Y, and Rk
%                W : basis for solution search space - dimension k = trunc_options.nOuter;
%                Y : basis for rhs search space
%               Rk : Rk-factor, AW = YR
%                    AW = YRk
%                    approx solution tilde(x) = Ws_k
%                x : initial solution for the case that the bases are
%                    given. x can be [] or some appropriate initial base
%
% Outputs:
%      x_out : computed solution
%     output : structure with the following fields:
%           x_all - all x in each iteration  
%      iterations - stopping iteration (options.Iter | GCV-determined)
%         GCVstop - GCV curve used to find stopping iteration
%            Enrm - relative error norms (requires x_true)
%            Rnrm - relative residual norms
%            Xnrm - relative solution norms
%            U,QV - U and V are genGK basis vectors and Q is from the prior
%               B - bidiagonal matrix from genGK
%            flag - a flag that describes the output/stopping condition:
%                       1 - flat GCV curve
%                       2 - min of GCV curve (within window of MinTol its)
%                       3 - performed max number of iterations
%                       4 - achieved residual tolerance
%           alpha - regularization parameter at (output.iterations) its
%           Alpha - vector of all regularization parameters computed
%     truncmats : structure containing recycled matrices, W, Y, Rk, and x
%
% References:
%   [1] Chung, Nagy and O'Leary, "A Weighted-GCV Method for Lanczos-Hybrid
%       Regularization", ETNA 28 (2008), pp. 149-167.
%   [2]  Chung, de Sturler, and Jiang. "Hybrid Projection Methods with
%           Recycling for Inverse Problems". SISC, 2020.
%
% J. Chung, E. de Sturler, J. Jiang, C. Ma, 2021

%% Initialization
defaultopt = struct('InSolv','tikhonov','RegPar','wgcv','Omega',...
  'adapt', 'Iter', [], 'Reorth','off','x_true', 'off', 'BegReg', 2,...
  'Vx' , [], 'FlatTol', 10^-6, 'MinTol', 4, 'ResTol', [10^-6, 10^-6]);

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
  x_out = defaultopt;
  return;
end

% Check for acceptable number of input arguments
if nargin < 4
  error('HyBR: Not Enough Inputs')
elseif nargin < 5
  P = {[],[]}; options = []; trunc_options = [];
elseif nargin < 6
  options = []; trunc_options = [];
end
if isempty(options)
  options = defaultopt;
end

if isempty(trunc_options) % Additional parameters for truncation
  nOuter = 5;
  nInner = 40;
else
  nOuter = trunc_options.nOuter;
  nInner = trunc_options.nInner;
end

if iscell(P) % Preconditioner
  P_le = P{1};
  P_ri = P{2};
else
  P_le = P;
  P_ri = [];
end

% Get options:
[m,n] = size(A);
defaultopt.Iter = min([m, n, 100]);
options = HyBRset(defaultopt, options);

solver = HyBRget(options,'InSolv',[],'fast');
regpar = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');
maxiter = HyBRget(options,'Iter',[],'fast');
x_true = HyBRget(options,'x_true',[],'fast');
regstart = HyBRget(options,'BegReg',[],'fast');
degflat = HyBRget(options,'FlatTol',[],'fast');
mintol = HyBRget(options,'MinTol',[],'fast');
restol = HyBRget(options,'ResTol',[],'fast');

adaptWGCV = strcmp(regpar, {'wgcv'}) && strcmp(omega, {'adapt'});
notrue = strcmp(x_true,{'off'});

%--------------------------------------------
%  The following is needed for RestoreTools:
%
if isa(A, 'psfMatrix')
  bSize = size(b);
  b = b(:);
  A.imsize = bSize;
  if ~notrue
    
    x_true = x_true(:);
  end
end
%
%  End of new stuff needed for RestoreTools
%--------------------------------------------

% Set-up output parameters:
outputparams = nargout>1;
if outputparams
  output.iterations = maxiter;
  output.GCVstop = [];
  output.Enrm = ones(maxiter,1);
  output.Rnrm = ones(maxiter,1);
  output.Xnrm = ones(maxiter,1);
  output.U = [];
  output.Q = [];
  output.QV = [];
  output.B = [];
  output.x_all = {};
  output.flag = 3;
  output.alpha = 0;
  output.E_opt = ones(maxiter,1);
end

% Test for a left preconditioner and define solver: 
if isempty(P_le)
  beta = normM(b,@(x)R\x);
  U = (1 / beta)*b;
  handle = @genGKB;
  if ~isempty(P_ri)
    handle = @genGKBP;     %preconditioner additional function?
  end
else
  U = P_le\b;
  beta = normM(U,@(x)R\x); U = U / beta;
  handle = @genGKBP;
end

switch solver
  case 'tsvd'
    solverhandle = @TSVDsolver;
  case 'tikhonov'
    solverhandle = @Tikhonovsolver;
end

% check/setup reycling
if isempty(trunc_mats)
  no_recycle = 1; % build recycling matrices at start
else
  no_recycle = 0;
  % W and S must be defined; check size
  if isempty(trunc_mats.W) | size(trunc_mats.W,2) > (nInner-1)
    fprintf('trunc_mats.W must exist \n');
    fprintf('number of recycling vectors must be less than max_mm+2 \n\n');
    return
  else
    W = trunc_mats.W; % W is assumed to have orthonormal columns
    kk = size(W,2);
    if isempty(trunc_mats.x) % no initial solution is given
      x_out = zeros(size(A,2),1);
    else                      % initial solution is given            (construction of W_k)
      x_out = trunc_mats.x;
      eta = W'*x_out; xhat = x_out - W*eta; etahat = norm(xhat);
      xhat = xhat / etahat;
      W = [W xhat];
      kk = kk+1;
    end
  end
  if isempty(trunc_mats.Y)                   % Skinny QR factorization, Y is orthonomal matrix and Rk is upper-triangular matrix
    Y = zeros(size(A,1),size(W,2));
    for j = 1:kk
      Y(:,j) = A*(Q*W(:,j));
    end
    [Y,Rk] = qr(Y,0);
  else
    Y = trunc_mats.Y;
    Rk = trunc_mats.Rk;
  end
end

%% Main Code Begins Here
B = []; V = []; QV = []; GCV = []; Omega = []; Alpha = [];
terminate = 0;
if (strcmp(regpar, {'wgcv'}) || strcmp(regpar, {'gcv'}))
  terminate = 1;
end
warning = 0; iter = 0;
  h = waitbar(0,'Performing genHyBRrecycle, please wait...');
for outer = 1:nOuter

  if no_recycle
    % first run is genGK with nInner iterations
    for inner = 1:nInner
      iter = iter+1;
      [U, B, V, QV] = feval(handle, A, Q, R, U, B, V, QV, P_le, P_ri, options);
      if strcmp(regpar,'optimal')
        options.Vx = V;
      end
      vector = (beta*eye(size(B,2)+1,1)); % assumes b is first vector in V
      if inner >= 2 % Otherwise skip
          switch solver
          case{'tsvd', 'tikhonov'}% Solve projected problem using TSVD/Tikhonov at each iteration
            [Ub, Sb, Vb] = svd(B);
            
            if adaptWGCV %Use the adaptive, weighted GCV method
              Omega(iter-1) = min(1, findomega(Ub'*vector, diag(Sb), solver));
              options.Omega = mean(Omega);
            end
            
            [y, alpha] = feval(solverhandle, Ub, diag(Sb), Vb, vector, options, B, beta, QV, m);
            Alpha(iter-1) = alpha;
            % Compute the GCV value used to find the stopping criteria
            GCV(iter-1) = GCVstopfun(alpha, Ub(1,:)', diag(Sb), beta, m, n, solver);
            
            % Determine if GCV wants us to stop
            if iter > 2 && terminate
              %%-------- If GCV curve is flat, we stop -----------------------
              if abs((GCV(iter-1)-GCV(iter-2)))/GCV(regstart-1) < degflat
                x_out = QV*y; % Return the solution at (i-1)st iteration
                % Test for a right preconditioner:
                if ~isempty(P_ri)
                  x_out = P_ri \ x_out;
                end
                if notrue %Set all the output parameters and return
                  if outputparams
                    output.U = U;
                    output.QV = QV;
                    output.V = V;
                    output.B = B;
                    output.GCVstop = GCV(:);
                    output.iterations = iter-1;
                    output.flag = 1;
                    output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                    output.Alpha = Alpha(1:iter-1); % Reg Parameters
                  end
                  close(h)
                  %--------------------------------------------
                  %  The following is needed for RestoreTools:
                  %
                  if isa(A, 'psfMatrix')
                    x_out = reshape(x_out, bSize);
                  end
                  %
                  %  End of new stuff needed for RestoreTools
                  %--------------------------------------------
                  return;
                else % Flat GCV curve means stop, but continue since have x_true
                  if outputparams
                    output.iterations = iter-1; % GCV says stop at (i-1)st iteration
                    output.flag = 1;
                    output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                  end
                end
                terminate = 0; % Solution is already found!
                
                %%--- Have warning : Avoid bumps in the GCV curve by using a
                %    window of (mintol+1) iterations --------------------
              elseif warning && length(GCV) > iterations_save + mintol %Passed window
                if GCV(iterations_save) < GCV(iterations_save+1:end)
                  % We should have stopped at iterations_save.
                  x_out = x_save;
                  % Test for a right preconditioner:
                  if ~isempty(P_ri)
                    x_out = P_ri \ x_out;
                  end
                  if notrue %Set all the output parameters and return
                    if outputparams
                      output.U = U;
                      output.V = V;
                      output.QV = QV;
                      output.B = B;
                      output.GCVstop = GCV(:);
                      output.iterations = iterations_save;
                      output.flag = 2;
                      output.alpha = alpha_save;
                      output.Alpha = Alpha(1:iterations_save); % Reg Parameters
                    end
                    close(h)
                    %--------------------------------------------
                    %  The following is needed for RestoreTools:
                    %
                    if isa(A, 'psfMatrix')
                      x_out = reshape(x_out, bSize);
                    end
                    %
                    %  End of new stuff needed for RestoreTools
                    %--------------------------------------------
                    return;
                  else % GCV says stop at iterations_save, but continue since have x_true
                    if outputparams
                      output.iterations = iterations_save;
                      output.flag = 2;
                      output.alpha = alpha_save;
                    end
                  end
                  terminate = 0; % Solution is already found!
                  
                else % It was just a bump... keep going
                  warning = 0;
                  x_out = [];
                  iterations_save = maxiter;
                  alpha_save = 0;
                end
                
                %% ----- No warning yet: Check GCV function---------------------
              elseif ~warning
                if GCV(iter-2) < GCV(iter-1) %Potential minimum reached.
                  warning = 1;
                  % Save data just in case.
                  x_save = QV*y;
                  iterations_save = iter-1;
                  alpha_save = alpha;
                end
              end
            end
            
          case 'none'
            y = B \ vector;
            alpha = 0;
          otherwise
            error('genHyBR error: No inner solver!')
        end
        ri = vector - B*y;
        ri_nrm = norm(ri);
        Vy = QV*y;
        r_out = b - A*Vy;
        rnrm = norm(r_out);
        
        if outputparams
          if ~notrue
            output.Enrm(iter-1,1) = norm(Vy(:)-x_true(:))/norm(x_true(:));
            x_opt = QV*(QV\x_true);                                 
            output.E_opt(iter-1,1) = norm( x_opt(:)-x_true(:) )/norm(x_true(:));
          end
          output.Rnrm(iter-1,1) = rnrm;
          output.Xnrm(iter-1,1) = norm(Vy);
          output.x_all{iter-1,1} = Vy;
        end
        
        if ri_nrm < restol(1)*beta
          Vy = QV*y; % FIX - already above
          x_out = Vy; % FIX - already above
          return
        end
      end
    end
    
    output.B = B;
    output.QV = QV;
    output.V = V;
    output.U = U;
    
    % --------- Select a recycle space W --------------------------
    y_1 = y; % regularized solution y
    Vy = QV*y_1;
    %%%%%%%%%%%%%% Perform Compression %%%%%%%%%%%%%%%%%%%
    [kk,W,Y,Rk] = compression(A,B,V,trunc_options,Vy,vector,y,[]);
    
    x_out = Vy;
    no_recycle = 0;
    
  else % no_reycle = false, but possible first outer iteration
    if outer == 1 % set up
      r_out = b - A*x_out;
    else
      
      Vupd = W(:,end);
      r_out = b - A*Vupd;
    end
    
    % build up space again
    % initialize next inner
    r_upd = r_out;
    
    zeta = (Y'*r_upd);
    btil = r_upd - Y*zeta;
    betaInn = normM(btil,@(x)R\x); U = btil/betaInn;
    B = []; V = []; H = []; GCV = []; Omega = []; QV = [];
    handle = @genrecyclingGKB; %% recyclingGKB process
    
    for inner = 1:nInner-kk
      iter = iter+1;
      [U, B, V, QV, H] = feval(handle, A, U, B, V, QV, Q, R, H, P_le, P_ri, W, Y, options);
      if strcmp(regpar,'optimal')
        options.Vx = [W,V];
      end
      mm = inner;
      BB = zeros(kk+mm+1,kk+mm); % we could extend BB rather then compute from scratch
      BB(1:kk,1:kk) = Rk;
      BB(1:kk,kk+1:kk+mm) = H;
      BB(kk+1:kk+mm+1,kk+1:kk+mm) = B;
      vector = zeros(mm+kk+1,1);
      vector(1:kk) = zeta + Rk(:,end);
      vector(kk+1:kk+mm+1) = (betaInn*eye(size(B,2)+1,1)); % assumes normalized btil is first vector in U
      
      [Ub,Sb,Vb] = svd(BB);
      
      switch solver
        case{'tsvd', 'tikhonov'}% Solve projected problem using TSVD/Tikhonov at each iteration
          if mm+kk == 1
            Sb = Sb(1,1);
            
            if adaptWGCV %Use the adaptive, weighted GCV method
              if iter > 1
                Omega(iter-1) = min(1, findomega(Ub'*vector,Sb, solver));
              else
                Omega(1) = 1;
              end
              options.Omega = mean(Omega);
            end
            
            y = feval(solverhandle, Ub, Sb, Vb, vector, options, BB, iter, [W,V],m);
                     
          else
            
            if adaptWGCV %Use the adaptive, weighted GCV method
              if iter>1
                Omega(iter-1) = min(1, findomega(Ub'*vector,diag(Sb), solver));
              else
                Omega(1) = 1;
              end
              options.Omega = mean(Omega);
            end
            [y, alpha] = feval(solverhandle, Ub, diag(Sb), Vb, vector, options, BB, iter, [W,V],m);
            
          end
          output.Alpha(iter) = alpha;
          
          % Compute the GCV value used to find the stopping criteria
          if iter > 1
            GCV(iter-1) = GCVstopfun(alpha, Ub(1,:)', diag(Sb), beta, m, n, solver);
          end
          % Determine if GCV wants us to stop
          if iter > 2 && terminate
            %%-------- If GCV curve is flat, we stop -----------------------
            if abs((GCV(iter-1)-GCV(iter-2)))/GCV(regstart-1) < degflat
              x_out = [W,V]*y; % Return the solution at (i-1)st iteration
              % Test for a right preconditioner:
              if ~isempty(P_ri)
                x_out = P_ri \ x_out;
              end
              if notrue %Set all the output parameters and return
                if outputparams
                  output.U = [Y,U];
                  output.V = [W,V];
                  output.B = BB;
                  output.x_all{iter-1,1} = [W,V]*y;
                  output.GCVstop = GCV(:);
                  output.iterations = iter-1;
                  output.flag = 1;
                  output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                  output.Alpha = Alpha(1:iter-1); % Reg Parameters
                end
                close(h)
                %--------------------------------------------
                %  The following is needed for RestoreTools:
                %
                if isa(A, 'psfMatrix')
                  x_out = reshape(x_out, bSize);
                end
                %
                %  End of new stuff needed for RestoreTools
                %--------------------------------------------
                return;
              else % Flat GCV curve means stop, but continue since have x_true
                if outputparams
                  output.iterations = iter-1; % GCV says stop at (i-1)st iteration
                  output.flag = 1;
                  output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                end
              end
              terminate = 0; % Solution is already found!
              
              %%--- Have warning : Avoid bumps in the GCV curve by using a
              %    window of (mintol+1) iterations --------------------
            elseif warning && length(GCV) > iterations_save + mintol %Passed window
              if GCV(iterations_save) < GCV(iterations_save+1:end)
                % We should have stopped at iterations_save.
                x_out = x_save;
                % Test for a right preconditioner:
                if ~isempty(P_ri)
                  x_out = P_ri \ x_out;
                end
                if notrue %Set all the output parameters and return
                  if outputparams
                    output.U = [Y,U];
                    output.V = [W,V];
                    output.B = BB;
                    output.GCVstop = GCV(:);
                    output.x_all{iter-1,1} = [W,V]*y;
                    output.iterations = iterations_save;
                    output.flag = 2;
                    output.alpha = alpha_save;
                    output.Alpha = Alpha(1:iterations_save); % Reg Parameters
                  end
                  close(h)
                  %--------------------------------------------
                  %  The following is needed for RestoreTools:
                  %
                  if isa(A, 'psfMatrix')
                    x_out = reshape(x_out, bSize);
                  end
                  %
                  %  End of new stuff needed for RestoreTools
                  %--------------------------------------------
                  return;
                else % GCV says stop at iterations_save, but continue since have x_true
                  if outputparams
                    output.iterations = iterations_save;
                    output.flag = 2;
                    output.alpha = alpha_save;
                  end
                end
                terminate = 0; % Solution is already found!
                
              else % It was just a bump... keep going
                warning = 0;
                iterations_save = maxiter;
                alpha_save = 0;
              end
              
              %% ----- No warning yet: Check GCV function---------------------
            elseif ~warning
              if GCV(iter-2) < GCV(iter-1) %Potential minimum reached.
                warning = 1;
                % Save data just in case.
                x_save = [W,V]*y;
                iterations_save = iter-1;
                alpha_save = alpha;
              end
            end
          end
        case 'none' % Solve projected problem with no regularization
          y = BB \ vector;
      end
      
      x1 = y(1:kk); x2 = y(kk+1:mm+kk); % x1 = c and x2 = d in notes
      ri = vector - BB*y;
      ri_nrm = norm(ri);
      
      if outer > 1
        c = x1;
        d = x2;
      else
        c = x1;
        d = x2;
      end
      
      x_out = W*c + V*d; % x_out is the adding vector
      r_out = b - A*x_out;
      rnrm = norm(r_out);
      
      if outputparams && iter > 1
        if ~notrue
          output.Enrm(iter-1,1) = norm(x_out(:)-x_true(:))/norm(x_true(:));
          % projection of x_true onto solution space
          x_opt = W*(W'*x_true) + V*(V'*x_true);
          output.E_opt(iter-1,1) = norm( x_opt(:) - x_true(:) )/norm(x_true(:));
        end
        output.Rnrm(iter-1,1) = rnrm;
        output.x_all{iter-1,1} = x_out;
        output.Xnrm(iter-1,1) = norm(x_out);
      end
      
      if ri_nrm < restol(1)*beta
        trunc_mats.W = W;
        trunc_mats.Rk = Rk;
        trunc_mats.Y = Y;
        output.iterations = iter;
        return
      end
      
    end
    
    % --------- Select a recycle space W (do compression) ---------------
    y_1 = [c;d]; % regularized solution, c is coefficient of input bases W, d is the coefficient of augmented bases
    
    y = y_1;
    x1 = y_1(1:kk); x2 = y_1(kk+1:mm+kk); 
    
    if outer > 1
      c = x1;
      d = x2;
    else
      c = x1;
      d = x2;
    end
    x_out = W*c + V*d;
    WV = [W V];
    %%%%%%%%%%%%%% Perform Compression %%%%%%%%%%%%%%%%%%%
    [kk,W,Y,Rk] = compression(A,BB,WV,trunc_options,x_out,vector,y,W);
  end
  waitbar(outer/nOuter,h);
end


close(h);
trunc_mats.W = W;
trunc_mats.Rk = Rk;
trunc_mats.Y = Y;
output.iterations = iter;
output.Enrm = output.Enrm(1:iter-1);
output.E_opt = output.E_opt(1:iter-1);
output.Rnrm = output.Rnrm(1:iter-1);
output.Xnrm = output.Xnrm(1:iter-1);


%% -----------------------SUBFUNCTION---------------------------------------
function omega = findomega(bhat, s, insolv)
%
%  
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the
%  derivative equal to zero and then solve for omega.
%
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%         insolv -  inner solver method for HyBR
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);
switch insolv
  case 'tsvd'
    k_opt = n;
    omega = (m*bhat(k_opt)^2) / (k_opt*bhat(k_opt)^2 + 2*bhat(k_opt+1)^2);
    
  case 'tikhonov'
    t0 = sum(abs(bhat(n+1:m)).^2);
    alpha = s(end);
    s2 = abs(s) .^ 2;
    alpha2 = alpha^2;
    
    tt = 1 ./ (s2 + alpha2);
    
    t1 = sum(s2 .* tt);
    t2 = abs(bhat(1:n).*alpha.*s) .^2;
    t3 = sum(t2 .* abs((tt.^3)));
    
    t4 = sum((s.*tt) .^2);
    t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);
    
    v1 = abs(bhat(1:n).*s).^2;
    v2 = sum(v1.* abs((tt.^3)));
    
    %
    % Now compute omega.
    %
    omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));
    
  otherwise
    error('Unknown solver');
end

%% ---------------SUBFUNCTION ---------------------------------------
function G = GCVstopfun(alpha, u, s, beta, m, n, insolv)
%
%  G = GCVstopfun(alpha, u, s, beta, n, insolv)
%  This function evaluates the GCV function G(i, alpha), that will be used
%     to determine a stopping iteration.
%
% Input:
%   alpha - regularization parameter at the kth iteration of HyBR
%       u - P_k^T e_1 where P_k contains the left singular vectors of B_k
%       s - singular values of bidiagonal matrix B_k
%    beta - norm of rhs b
%     m,n - size of the ORIGINAL problem (matrix A)
%  insolv - solver for the projected problem
%

k = length(s);
beta2 = beta^2;

switch insolv
  case 'tsvd'
    t2 = (abs(u(alpha+1:k+1))).^2;
    G = n*beta2*(sum(t2))/((m - alpha)^2);
  case 'tikhonov'
    s2 = abs(s) .^ 2;
    alpha2 = alpha^2;
    
    t1 = 1 ./ (s2 + alpha2);
    t2 = abs(alpha2*u(1:k) .* t1) .^2;
    t3 = s2 .* t1;
    
    num = beta2*(sum(t2) + abs(u(k+1))^2)/n;
    den = ( (m - sum(t3))/n )^2;
    G = num / den;
    
  otherwise
    error('Unknown solver');
end

