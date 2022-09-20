function [P,N] = All_HyBRs(LABEL,A,b,Q1,Q2,R,HYBR_TYPE,INSOLVE,REG,X_TRUE,nLevel,MAXIT,P,N,mu,trunc_options,trunc_mats)
%
% [P,N] = All_HyBRs(LABEL,A,b,Q1,Q2,R,HYBR_TYPE,INSOLVE,REG,X_TRUE,nLevel,MAXIT,P,N,mu)
%
% This function is to run and compare differen type of hybrid methods 
% combining subspace and regularization methods.
%
%   Input:
%       LABEL : Give a label to inidicate in Plot
%           A : either (a) a full or sparse matrix
%                     (b) a matrix object that performs mat*vec and mat'*vec
%              Note: If A is a function handle, create an object called funMat
%              e.g., A = funMat(Afun, Atfun) where Afun and Atfun are function
%                                            handles for A*vec and A'*vec
%           b : rhs vector
%       Q1,Q2 : covariance matrices, Q1 and Q2 is either (a) a full or sparse matrix or
%               (b) a matrix object that performs matrix*vector operations
%               At least Q1 has to be symmetric positive definite.
%               Q2 can be symmetric positive semidefinite.
%           R : covariance matrices of noise in data. (diagonal so easy to invert)
%   HYBR_TYPE : Hybrid methods [ LSQR | HyBR | genHyBR | mixHyBR ]          
%               or create yourself
%     INSOLVE : solver for the inner problem: [ none | TSVD | {Tikhonov} ] 
%               INSOLVE is only for LSQR, HyBR, and genHyBR
%         REG : solver for the inner problem: [ optimal | upre | gcv | wgcv ]
%               Note: 'optimal' requires x_true
%      X_TRUE : True solution : [ array | {off} ]
%               Returns error norms with respect to x_true at each iteration
%               and is used to compute 'optimal' regularization parameters
%      nLevel : if InSolv is 'upre', then nLevel represents the noise level 
%               and must be [non-negative scalar | {est}]
%         mu : prior mean of x
%              (default is zero vector. The (sample) prior mean is recommended)
%           P : P is a structure storing all computed results from the selected
%               hybrid method.
%           N : N is an index of structure to store in the structure P.
%
%   Outputs:
%
%           P : P is a structure storing all computed results from the selected
%               hybrid method.
%           N : N is an index of structure to store in the structure P.
%
% T. Cho 10/2019


switch HYBR_TYPE

    case 'HyBR'
        if isempty(nLevel)
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        else
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT);
        end
        [x, output] = HyBR(A, b(:), [], input);
        
    case 'HyBR0'
        if isempty(nLevel)
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        else
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT);
        end
        [x, output] = HyBR_0(A, b(:), mu, [], input);
        
    case 'LSQR'
        input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        [x, output] = HyBR(A, b(:), [], input);
    case 'genHyBR'
        if isempty(nLevel)
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        else
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT);
        end
        [x, output] = genHyBR(A, b(:), Q1, R, input);
        
    case 'genHyBR0'
        if isempty(nLevel)
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        else
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT);    
        end
        [x, output] = genHyBR_0(A, b(:), Q1, R,mu, input);
        
    case 'mixHyBR'
        if isempty(nLevel)
            input = mixHyBRset('RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT,'ReOrth', 'on','ResTol', [10^-6, 10^-6]);
        else
            input = mixHyBRset('RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT,'ReOrth', 'on','ResTol', [10^-6, 10^-6]);
        end
        [x, output] = mixHyBR(A, b(:), Q1, Q2, R, mu, input);
    
    case 'mixHyBRapprox0'
        if isempty(nLevel)
            input = mixHyBRset('RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT,'ReOrth', 'on');
        else
            input = mixHyBRset('RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT,'ReOrth', 'on');
        end
        [x, output] = mixHyBR_approx0(A, b(:), Q1, Q2, R,mu,input);
    case 'genHyBRrecycle'
        if isempty(nLevel)
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'Iter', MAXIT);
        else
            input = HyBRset('InSolv',INSOLVE,'RegPar',REG, 'x_true', X_TRUE(:),'nLevel',nLevel,'Iter', MAXIT);
        end
        [x, output,trunc_mats] = genHyBRrecycle(A, b(:), [], Q1, R, input, trunc_options, trunc_mats);
        
        
end

N = N + 1;
P{N,1} = LABEL;
% P{N,2} = HYBR_TYPE; 
P{N,2} = REG; 
P{N,3} = x; 
P{N,4} = output;
P{N,5} = trunc_mats;
    
    
