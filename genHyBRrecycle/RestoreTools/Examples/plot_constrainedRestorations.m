%
%  SCRIPT: plot_constrainedRestorations
%
%  Plot the restoration results of WMRNSD, KWMRNSD, and unconstrained LSQR
%  for the test problems used in IR_ConstrainedExamples script

subplot(1,3,1)
    imshow(x_lsqr, [0, max(x_lsqr(:))])
    title('LSQR')
subplot(1,3,2)
    imshow(x_wmrnsd, [0, max(x_wmrnsd(:))])
    title('WMRNSD')
subplot(1,3,3)
    imshow(x_kwmrnsd, [0, max(x_kwmrnsd(:))])
    title('KWMRNSD')