%
%  SCRIPT: IR_ConstrainedExamples
%
%  These examples show the basic usage of constrained iterative methods
%  (described in [1]) on some sample test data from RestoreTools. In
%  particular these codes reproduce the numerical results presented in
%  section 6.2 of [1].
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
%  The first example is the convergence history of the nonnegatively
%  constrained iterative methods for the"AtmosphericBlur30" test problem from 
%  RestoreTools [2].
disp(' ')
disp('Press any key to see the convergence history of the nonnegatively')
disp('constrained iterative methods for the AtmosphericBlur30 test problem')
pause
%
%   Load data
load AtmosphericBlur30
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'zero');
%
%  Now compare a few nonnegatively constrained iterative methods. We'll use 
%  mostly default options, but if we want to look at relative errors using 
%  the truth, and run the Richardson iteration using a default choice for 
%  tau (tau = 1/(||K||_1 ||K||_inf)) we need to first use IRset as follows:
options = IRset('x_true', f_true,'StepLength',-1);
%
%   Use the nonngatively constrained Richardson iteration to solve
[x_rinn, IterInfo_rinn] = IRgdnn(K, g, options);
[m_rinn, k_rinn] = min(IterInfo_rinn.Enrm);
%
%   Use nonnegatively constrained steepest descent to solve
options = IRset(options,'StepLength',0);
[x_gdnn, IterInfo_gdnn] = IRgdnn(K, g, options);
[m_gdnn, k_gdnn] = min(IterInfo_gdnn.Enrm);
%
%  Now let's try Modified Residual Norm Steepest Descent
[x_mrnsd, IterInfo_mrnsd] = IRmrnsd(K, g, options);
[m_mrnsd, k_mrnsd] = min(IterInfo_mrnsd.Enrm);
%
%  Compare with Weighted Modified Residual Norm Steepest Descent
[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_wmrnsd, k_wmrnsd] = min(IterInfo_wmrnsd.Enrm);
%
%  Compare with k-Weighted Modified Residual Norm Steepest Descent
[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_kwmrnsd, k_kwmrnsd] = min(IterInfo_kwmrnsd.Enrm);
%
%  Try Richardson-Lucy algorithm
[x_rl, IterInfo_rl]  = IRrl(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_rl, k_rl] = min(IterInfo_rl.Enrm);
%
%   Display convergence results
disp('(see Figure 1)')
figure(1), clf
axes('FontSize',FS)
plot( IterInfo_rinn.Enrm,'-ro','LineWidth',LW,'MarkerSize',MS), hold on
xlim([0 MaxIter+5])
plot( IterInfo_gdnn.Enrm,'-s','Color', myGreen,'LineWidth',LW,'MarkerSize',MS ), hold on
plot( IterInfo_mrnsd.Enrm,'-bd','LineWidth',LW,'MarkerSize',MS)
plot( IterInfo_wmrnsd.Enrm,'-m^','LineWidth',LW,'MarkerSize',MS)
plot( IterInfo_kwmrnsd.Enrm,'-cV','LineWidth',LW,'MarkerSize',MS)
plot( IterInfo_rl.Enrm,'-k>','LineWidth',LW,'MarkerSize',MS)
legend('RINN', 'SDNN', 'MRNSD', 'WMRNSD', 'KWMRNSD','RL','Location','Best');
delete(findobj(gca,'Type','line','Marker','o'));
delete(findobj(gca,'Type','line','Marker','s'));
delete(findobj(gca,'Type','line','Marker','d'));
delete(findobj(gca,'Type','line','Marker','^'));
delete(findobj(gca,'Type','line','Marker','V'));
delete(findobj(gca,'Type','line','Marker','>'));
plot( IterInfo_rinn.Enrm,'r','LineWidth',LW)
plot( IterInfo_gdnn.Enrm,'Color', myGreen,'LineWidth',LW )
plot( IterInfo_mrnsd.Enrm,'b','LineWidth',LW)
plot( IterInfo_wmrnsd.Enrm,'m','LineWidth',LW)
plot( IterInfo_kwmrnsd.Enrm,'c','LineWidth',LW)
plot( IterInfo_rl.Enrm,'k','LineWidth',LW)
plot(1:11:MaxIter,IterInfo_rinn.Enrm(1:11:MaxIter),'ro','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_gdnn.Enrm(1:11:MaxIter),'s','Color',myGreen,'MarkerSize',MS, 'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_mrnsd.Enrm(1:11:MaxIter),'bd','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_wmrnsd.Enrm(1:11:MaxIter),'m^','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_kwmrnsd.Enrm(1:11:MaxIter),'cV','MarkerSize',MS,'LineWidth',LW)
plot(1:11:MaxIter,IterInfo_rl.Enrm(1:11:MaxIter),'k>','MarkerSize',MS,'LineWidth',LW)
xlabel('Iterations','fontsize',FS)
ylabel('Relative error','fontsize',FS)


disp(' ')
disp('Press any key to see results of WMRNSD, KWMRNSD, and unconstrained')
disp('LSQR for the AtmosphericBlur30 test problem')
pause
%
%  Now let's try lsqr:
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);
%
%  Display convergence results
disp('(see Figure 2)')
figure(2), clf
plot_constrainedConvergence
%
%  Compute restorations at 47th iteration    
minIter = 47;
options = IRset(options,'MaxIter', minIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
%
%  Display restoration results
disp('(see Figure 3)') 
figure(3), clf
plot_constrainedRestorations


%------------------------------------------------------------------------
%
% The second example is the "GaussianBlur422"test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of WMRNSD, KWMRNSD, and unconstrained')
disp('LSQR for the GaussianBlur422 test problem')
pause
%
%   Load data
load GaussianBlur422
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'reflexive');
%
%  Now compare WMRNSD and KWMRNSD with the unconstrained LSQR
options = IRset(options, 'x_true', f_true,'MaxIter',MaxIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);

[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_wmrnsd, k_wmrnsd] = min(IterInfo_wmrnsd.Enrm);

[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_kwmrnsd, k_kwmrnsd] = min(IterInfo_kwmrnsd.Enrm);
%
%  Display convergence results
disp('(see Figure 4)')
figure(4), clf
plot_constrainedConvergence
%
%  Compute restorations at 29th iteration    
minIter = 29;
options = IRset(options,'MaxIter', minIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
%
%  Display restoration results
disp('(see Figure 5)') 
figure(5), clf
plot_constrainedRestorations
       
%------------------------------------------------------------------------
%
%  The third example is the "VariantGaussianBlur2" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of WMRNSD, KWMRNSD, and unconstrained')
disp('LSQR for the VariantGaussianBlur2 test problem')
pause
%
%  Load data
load VariantGaussianBlur2
%
%  For this data we need to construct a psfMatrix.  
K = psfMatrix(PSF, center,'zero');
%
%  Now compare WMRNSD and KWMRNSD with the unconstrained LSQR
options = IRset(options, 'x_true', f_true,'MaxIter',MaxIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);

[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_wmrnsd, k_wmrnsd] = min(IterInfo_wmrnsd.Enrm);

[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_kwmrnsd, k_kwmrnsd] = min(IterInfo_kwmrnsd.Enrm);
%
%  Display convergence results
disp('(see Figure 6)')
figure(6), clf
plot_constrainedConvergence

%
%  Compute restorations at 20th iteration
minIter = 20;
options = IRset(options,'MaxIter', minIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
%
%  Display restoration results
disp('(see Figure 7)') 
figure(7), clf
plot_constrainedRestorations
    
%------------------------------------------------------------------------
%
% The third example is the "VariantMotionBlur_medium" test problem from 
% RestoreTools [2].
disp(' ')
disp('Press any key to see results of WMRNSD, KWMRNSD, and unconstrained')
disp('LSQR for the VariantMotionBlur_medium test problem')
pause
%
%   Load data
load VariantMotionBlur_medium
%
%  Now compare WMRNSD and KWMRNSD with the unconstrained LSQR
options = IRset(options, 'x_true', f_true, 'MaxIter', MaxIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[m_lsqr, k_lsqr] = min(IterInfo_lsqr.Enrm);

[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_wmrnsd, k_wmrnsd] = min(IterInfo_wmrnsd.Enrm);

[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[m_kwmrnsd, k_kwmrnsd] = min(IterInfo_kwmrnsd.Enrm);

%
%  Display convergence results
disp('(see Figure 8)')
figure(8), clf
plot_constrainedConvergence
%
%  Compute restorations at 35th iteration
minIter = 35;
options = IRset(options,'MaxIter', minIter);
[x_lsqr, IterInfo_lsqr] = IRlsqr(K, g, options);
[x_wmrnsd, IterInfo_wmrnsd]  = IRwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
[x_kwmrnsd, IterInfo_kwmrnsd]  = IRkwmrnsd(K, g, options, NoiseInfo.GaussianStdev^2, NoiseInfo.PoissonParam);
%
%  Display restoration results
disp('(see Figure 9)') 
figure(9), clf
plot_constrainedRestorations
    

