%
%  SCRIPT: IR_UnconstrainedExamples
%
%  These examples show the basic usage of unconstrained iterative methods
%  (described in [1]) on some sample test data from RestoreTools. In
%  particular these codes reproduce the numerical results presented in
%  section 6.1 of [1].
%
%
%  References: 
%  [1] S. Berisha and J. Nagy, "Iterative Methods for Image Restoration"
%  [2] J. Nagy, K. Palmer, L. Perrone, "Iterative Methods for Image 
%      Deblurring: A Matlab Object Oriented Approach",
%      Numerical Algorithms, 36 (2004), pp. 73-93.
%      http://www.mathcs.emory.edu/~nagy/RestoreTools
%
%  S. Berisha and J. Nagy, February 2012
%

FS = 14;
LW = 1.5;
MS = 10;
MaxIter = 100;
myGreen = [0 0.5 0];


%------------------------------------------------------------------------
%
% The first example is the "AtmosphericBlur30" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of AtmosphericBlur30 test problem')
pause
%
%   Load data
load AtmosphericBlur30
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'zero');
%
%  Now compare a few basic unconstrained iterative methods. We'll use 
%  mostly default options, but if we want to look at relative errors using 
%  the truth, and run the Richardson iteration using a default choice for 
%  tau (tau = 1/(||K||_1 ||K||_inf)) we need to first use IRset as follows:
options = IRset('x_true', f_true,'StepLength',-1);
%
%  Use Richardson iteration to solve
[x_ri, IterInfo_ri] = IRgd(K, g, options);
[m_ri, k_ri] = min(IterInfo_ri.Enrm);
%
%  Use steepest descent to solve
options = IRset(options,'StepLength',0);
[x_sd, IterInfo_sd] = IRgd(K, g, options);
[m_sd, k_sd] = min(IterInfo_sd.Enrm);
%
%  Now let's try lsqr:
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);
%
%  Compare with cgls:
[x_cgls, IterInfo_cgls] = IRcgls(K, g, options);
[m_cgls, k_cgls] = min(IterInfo_cgls.Enrm);
%
%  Try Hybr
[x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
[m_hybr, k_hybr] = min(IterInfo_hybr.Enrm);
%
%  Display convergence results
disp('(see Figure 1)')
figure(1), clf
plot_unconstrainedConvergence
%
%  Calculate the overall minimum number of iterations it took for any of 
%  the methods to reach their minimum restoration error. Compute solutions 
%  corresponding to that number of iterations.
minIter = min([k_ri,k_sd, k_cgls, k_lsqr, k_hybr]) - 1;
if(minIter <100)
    options = IRset(options,'MaxIter', minIter,'StepLength',-1);
    [x_ri, IterInfo_ri] = IRgd(K, g, options);
    options = IRset(options,'MaxIter', minIter,'StepLength',0);
    [x_sd, IterInfo_sd] = IRgd(K, g, options);
    [x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
    [x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
end
%
%  Display restoration results
disp('(see Figure 2)') 
figure(2), clf
plot_unconstrainedRestorations
    

%------------------------------------------------------------------------
%
% The second example is the "GaussianBlur422"test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of GaussianBlur422 test problem')
pause
%
%   Load data
load GaussianBlur422
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'reflexive');
%
%  Now compare a few basic unconstrained iterative methods. We'll use 
%  mostly default options, but if we want to look at relative errors using 
%  the truth, and run the Richardson iteration using a default choice for 
%  tau (tau = 1/(||K||_1 ||K||_inf)) we need to first use IRset as follows:
options = IRset(options, 'x_true', f_true,'MaxIter', MaxIter, 'StepLength',-1 );
%
%   Use Richardson iteration to solve
[x_ri, IterInfo_ri] = IRgd(K, g, options);
[m_ri, k_ri] = min(IterInfo_ri.Enrm);
options = IRset(options,'StepLength',0);
%
%   Use steepest descent to solve
[x_sd, IterInfo_sd] = IRgd(K, g, options);
[m_sd, k_sd] = min(IterInfo_sd.Enrm);
%
%   Now let's try lsqr:
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);
%
%  Compare with cgls:
[x_cgls, IterInfo_cgls] = IRcgls(K, g, options);
[m_cgls, k_cgls] = min(IterInfo_cgls.Enrm);
%
%  Try Hybr
[x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
[m_hybr, k_hybr] = min(IterInfo_hybr.Enrm);

%
%   Display convergence results
disp('(see Figure 3)')
figure(3), clf
plot_unconstrainedConvergence
%
%  Calculate the overall minimum number of iterations it took for any of 
%  the methods to reach their minimum restoration error. Compute solutions 
%  corresponding to that number of iterations.
minIter = min([k_ri,k_sd, k_cgls, k_lsqr, k_hybr]) - 1;

if(minIter <100)
    options = IRset(options,'MaxIter', minIter,'StepLength',-1);
    [x_ri, IterInfo_ri] = IRgd(K, g, options);
    options = IRset(options,'MaxIter', minIter,'StepLength',0);
    [x_sd, IterInfo_sd] = IRgd(K, g, options);
    [x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
    [x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
end
%
%  Display restoration results
disp('(see Figure 4)') 
figure(4), clf
plot_unconstrainedRestorations
    
%------------------------------------------------------------------------
%
% The third example is the "VariantGaussianBlur2" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of VariantGaussianBlur2 test problem')
pause
%
%   Load data
load VariantGaussianBlur2
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'zero');
%
%  Now compare a few basic unconstrained iterative methods. We'll use 
%  mostly default options, but if we want to look at relative errors using 
%  the truth, and run the Richardson iteration using a default choice for 
%  tau (tau = 1/(||K||_1 ||K||_inf)) we need to first use IRset as follows:
options = IRset('x_true', f_true,'MaxIter', MaxIter, 'StepLength',-1);
%
%   Use Richardson iteration to solve
[x_ri, IterInfo_ri] = IRgd(K, g, options);
[m_ri, k_ri] = min(IterInfo_ri.Enrm);
options = IRset(options,'StepLength',0);
%
%   Use steepest descent to solve
[x_sd, IterInfo_sd] = IRgd(K, g, options);
[m_sd, k_sd] = min(IterInfo_sd.Enrm);
%
%   Now let's try lsqr:
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);
%
%  Compare with cgls:
[x_cgls, IterInfo_cgls] = IRcgls(K, g, options);
[m_cgls, k_cgls] = min(IterInfo_cgls.Enrm);
%
%  Try HyBR
[x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
[m_hybr, k_hybr] = min(IterInfo_hybr.Enrm);

%
%  Display convergence results
disp('(see Figure 5)')
figure(5), clf
plot_unconstrainedConvergence

minIter = min([k_ri,k_sd, k_cgls, k_lsqr, k_hybr]) - 1;

if(minIter <100)
    options = IRset(options,'MaxIter', minIter,'StepLength',-1);
    [x_ri, IterInfo_ri] = IRgd(K, g, options);
    options = IRset(options,'MaxIter', minIter,'StepLength',0);
    [x_sd, IterInfo_sd] = IRgd(K, g, options);
    [x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
    [x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
end
%
%  Display restoration results
disp('(see Figure 6)') 
figure(6), clf
plot_unconstrainedRestorations


%------------------------------------------------------------------------
%
% The fourth example is the "VariantMotionBlur_medium" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of VariantMotionBlur_medium test problem')
pause
%
%   Load data
load VariantMotionBlur_medium
%
%  Now compare a few basic unconstrained iterative methods. We'll use 
%  mostly default options, but if we want to look at relative errors using 
%  the truth, and run the Richardson iteration using a default choice for 
%  tau (tau = 1/(||K||_1 ||K||_inf)) we need to first use IRset as follows:
options = IRset('x_true', f_true, 'MaxIter', MaxIter, 'StepLength',-1);
%
%   Use Richardson iteration to solve
[x_ri, IterInfo_ri] = IRgd(K, g, options);
[m_ri, k_ri] = min(IterInfo_ri.Enrm);
options = IRset(options,'StepLength',0);
%
%   Use steepest descent to solve
[x_sd, IterInfo_sd] = IRgd(K, g, options);
[m_sd, k_sd] = min(IterInfo_sd.Enrm);
%
%   Now let's try lsqr:
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);
%
%  Compare with cgls:
[x_cgls, IterInfo_cgls] = IRcgls(K, g, options);
[m_cgls, k_cgls] = min(IterInfo_cgls.Enrm);
%
%  Try HyBR
[x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
[m_hybr, k_hybr] = min(IterInfo_hybr.Enrm);

%
%  Display convergence results
disp('(see Figure 7)')
figure(7), clf
plot_unconstrainedConvergence

minIter = min([k_ri,k_sd, k_cgls, k_lsqr, k_hybr]) - 1;

if(minIter <100)
    options = IRset(options,'MaxIter', minIter,'StepLength',-1);
    [x_ri, IterInfo_ri] = IRgd(K, g, options);
    options = IRset(options,'MaxIter', minIter,'StepLength',0);
    [x_sd, IterInfo_sd] = IRgd(K, g, options);
    [x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
    [x_hybr, IterInfo_hybr] = IRhybr(K, g, options);
end
%
%  Display restoration results
disp('(see Figure 8)') 
figure(8), clf
plot_unconstrainedRestorations
