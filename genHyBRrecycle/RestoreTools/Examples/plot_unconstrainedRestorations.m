%
%  SCRIPT: plot_unconstrainedRestorations
%
%  Plot the restoration results for the test problems used in
%  IR_UnconstrainedExamples script

subplot(2,2,1)
    imshow(x_ri, [0, max(x_ri(:))])
    title('RI')
subplot(2,2,2)
    imshow(x_sd, [0, max(x_sd(:))])
    title('SD')
subplot(2,2,3)
    imshow(x_lsqr, [0, max(x_lsqr(:))])
    title('LSQR')
subplot(2,2,4)
    imshow(x_hybr, [0, max(x_hybr(:))])  
    title('HyBR')
    